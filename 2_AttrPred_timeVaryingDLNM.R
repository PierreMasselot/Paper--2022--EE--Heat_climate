#################################################################
#
#           AUDACE - Summer mortality prediction
#         Year specific Heat attributable fraction
#
#################################################################

#----------------------------
# Apply DLNM and compute AF
#----------------------------

#----- Prepare objects to store results

# Store crosspred
cplist <- vector("list", length(region))
names(cplist) <- region

# Store MMT, confidence intervals
mmt <- matrix(NA, nrow = ny, ncol = length(region),
  dimnames = list(years, region))
mmtcis <- afcis <- meanafcis <- vector("list", length(region))
names(mmtcis) <- names(afcis) <- region

# Parallelization
ncores <- detectCores()
cl <- makeCluster(max(1, ncores - 2))
registerDoParallel(cl)

# Loop on regions
for (i in seq_along(dlist)){
  
  #----- Run time-varying DLNM
  # Data for region
  daily_dat <- dlist[[i]]
  
  # Compute temperature distribution percentiles for AF threshol
  regper <- quantile(daily_dat$Tmean, heatper, na.rm = T)
  
  # Specify crossbasis (see Gasparrini et al. 2015 EHP)
  cb <- crossbasis(daily_dat$Tmean, lag = maxlag, 
    argvar = list(fun = varfun, degree = vardeg, 
      knots = quantile(daily_dat$Tmean, varper)), 
    arglag = list(knots = logknots(maxlag, lagnk)))
  
  # Define day as numeric and scale
  datenum <- drop(scale(as.numeric(daily_dat$Date)))
  
  # Define interaction term
  int <- cb * datenum
  
  # Run the GLM
  res <- glm(Death ~ cb + dow + ns(doy, df = dfseas):factor(year) + 
      ns(Date, df = round(dfdec * ny / 10)) + int,
    data = daily_dat, family = "quasipoisson")
  
  #----- Extract crosspred at each first of July
  
  # Get scaled date for 1st of Julys
  july1 <- datenum[format(daily_dat$Date, "%m-%d") == "07-01"]
  
  # Indices to extract main coefs and interaction coef
  indmain <- 1:ncol(cb) + 1
  indint <- grep("int", names(coef(res)))
  
  # Extract coefs and combine for each year
  coefmain <- coef(res)[indmain]
  coefint <- coef(res)[indint]
  coefsets <- lapply(july1, function(x) coefmain + x * coefint)
  
  # Extract vcovs and combine for each year
  vcovmain <- vcov(res)[indmain, indmain]
  vcovint <- vcov(res)[indint, indint]
  vcovmainint <- vcov(res)[indmain, indint]
  vcovintmain <- vcov(res)[indint, indmain]
  vcovsets <- lapply(july1, function(x) vcovmain + x^2 * vcovint + 
      x * vcovintmain + x * vcovmainint)
  
  # Compute temperature distribution percentiles for crosspred
  tmeanper <- quantile(daily_dat$Tmean, predper / 100, na.rm = T)
  
  # Compute first crosspreds
  cps <- Map(crosspred, basis = list(cb), coef = coefsets,
    vcov = vcovsets, model.link = "log", at = list(tmeanper))
  
  # Determine MMT and update crosspred
  per <- quantile(daily_dat$Tmean, mmtper)
  inrange <- tmeanper >= per[1] & tmeanper <= per[2]
  mmt[,i] <- sapply(cps, function(x) 
    x$predvar[inrange][which.min(x$allRRfit[inrange])])
  cplist[[i]] <- Map(crosspred, basis = list(cb), coef = coefsets,
    vcov = vcovsets, model.link = "log", at = list(tmeanper), cen = mmt[,i])
  
  #----- Compute daily AN
  
  # Create centered basis for overall cumulative association 
  #   with centering changing each year
  bvar <- do.call(onebasis, c(list(x = daily_dat$Tmean), attr(cb, "argvar")))
  cenvec <- lapply(mmt[,i], function(cen) do.call(onebasis, 
    c(list(x = cen), attr(cb, "argvar"))))
  bvarcen <- Map(scale, x = split(as.data.frame(bvar), daily_dat$year),
    center = cenvec, scale = F)
  
  # MA of deaths for forward AN
  deaths <- lapply(split(daily_dat$Death, daily_dat$year), function(d)
    rowMeans(as.matrix(Lag(d, -maxlag:0)))
  )
  
  # Coefs for overall cumulative association
  creds <- Map(crossreduce, basis = list(cb), coef = coefsets,
    vcov = vcovsets, model.link = "log", at = list(tmeanper), cen = mmt[,i])
  totcoef <- lapply(creds, coef)
  
  # Compute contribution of each day
  anfun <- function(basis, coefs, deaths) 
    pmax((1 - exp(-basis %*% coefs)) * deaths, 0)
  anday <- Map(anfun, bvarcen, totcoef, deaths)
  
  #----- Compute AF
  
  # Above MMT
  datatab[[sprintf("AF_mmt_%s",  region[i])]] <- mapply(
    function(ydat, anday, mmt){
      sum(anday[ydat$Tmean > mmt], na.rm = T) / sum(ydat$Death)
    }, split(daily_dat, daily_dat$year), anday, mmt[,i])
  
  # Above temperature percentiles
  for (j in seq_along(heatper)){
    datatab[[sprintf("AF_%s_%s", heatper[j] * 100, region[i])]] <- 
      mapply(function(ydat, anday){
          sum(anday[ydat$Tmean > regper[j]], na.rm = T) / sum(ydat$Death)
        }, split(daily_dat, daily_dat$year), anday)
  }
  
  #----- Monte Carlo simulations for eCI
  
  # Generate new coefficients for GLM fit
  set.seed(1234 + i)
  coefregsimu <- mvrnorm(ns, coef(res), vcov(res))
  
  # Loop on simulations
  simures <- foreach(k = seq_len(ns), .packages = c("dlnm"), 
    .combine = abind) %dopar% {
    
    # Compute annual coefs
    coefk <- lapply(july1, function(x) coefregsimu[k, indmain] + 
        x * coefregsimu[k, indint])
    
    # Extract MMTs
    cps <- Map(crosspred, basis = list(cb), coef = coefk,
      vcov = vcovsets, model.link = "log", at = list(tmeanper))
    mmtk <- sapply(cps, function(x) 
      x$predvar[inrange][which.min(x$allRRfit[inrange])])
    
    # Reduce coefs
    creds <- Map(crossreduce, basis = list(cb), coef = coefk,
      vcov = vcovsets, model.link = "log", at = list(tmeanper), cen = mmtk)
    totcoef <- lapply(creds, coef)
    
    # Create centred basis
    cenvec <- lapply(mmtk, function(cen) do.call(onebasis, 
      c(list(x = cen), attr(cb, "argvar"))))
    bvarcen <- Map(scale, x = split(as.data.frame(bvar), daily_dat$year),
      center = cenvec, scale = F)
    
    # Compute contribution of each day
    anfun <- function(basis, coefs, deaths) 
      pmax((1 - exp(-basis %*% coefs)) * deaths, 0)
    anday <- Map(anfun, bvarcen, totcoef, deaths)
    
    # Compute AF
    afmmt <- mapply(
      function(ydat, anday, mmt){
        sum(anday[ydat$Tmean > mmt], na.rm = T) / sum(ydat$Death)
      }, split(daily_dat, daily_dat$year), anday, mmtk)
    afper <- sapply(regper, function(x) mapply(function(ydat, anday){
      sum(anday[ydat$Tmean > x], na.rm = T) / sum(ydat$Death)
    }, split(daily_dat, daily_dat$year), anday))
    
    # Output
    out <- cbind(mmtk, afmmt, afper)
    dim(out) <- c(dim(out), 1)
    out
  }
  
  #----- Compute confidence intervals
  
  # Compute quantiles of the array
  allcis <- apply(simures, 1:2, quantile, c(.025, .975))
  
  # Extract cis
  mmtcis[[i]] <- t(allcis[,,1])
  colnames(mmtcis[[i]]) <- c("low", "high")
  afcis[[i]] <- allcis[,,-1]
  
  # CI of average AF
  meanafcis[[i]] <- apply(apply(simures[,-1,], 2:3, mean), 1, 
    quantile, c(.025, .975))
}
  
# Stop parallel
stopCluster(cl)