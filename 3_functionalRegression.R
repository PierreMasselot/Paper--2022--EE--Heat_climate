#################################################################
#
#           AUDACE - Summer mortality prediction
#              Functional regression of AF
#
#################################################################

# Prepare parallelization
ncores <- detectCores()
cl <- makeCluster(max(1, ncores - 2))
registerDoParallel(cl)

#----------------------------
# Predict with climate indices for heat-related mortality defined through MMT
#----------------------------

#----- Loop on regions

# Loop on regions
reslist <- foreach(i = seq_along(region), .packages = c("FDboost")) %:%
  
  # Loop on outcomes
  foreach(j = seq_len(length(heatper) + 1), .packages = c("FDboost")) %dopar% {
    
    # Initialize list of results
    res <- vector("list", 2)
    names(res) <- c("main", "temp")
    
    #----- Main model
    # Prepare formula
    main_form <- sprintf("AF_%s_%s ~ %s + bbs(Year)", 
      c("mmt", heatper * 100)[j], region[i], 
      paste(sprintf("bsignal(%s, s, boundary.knots = c(%i, %i))", 
        climate_labs, clag[1] - 1, tail(clag, 1)), collapse = " + "))
    
    # Fit the model
    fdres <- FDboost(as.formula(main_form), timeformula = NULL, data = datatab,
      family = Gaussian(), control = boost_control(nu = nu, mstop = maxstop))
    
    # Select number of boosting steps
    cvres <- cvrisk(fdres, grid = 10:500,
      folds = cv(model.weights(fdres), type = "kfold", B = 10))
    mres <- mstop(cvres)
    res$main$cv <- cvres
    res$main$mstop <- mres
    
    # Selection of best model and kept variables
    res$main$final_mod <- fdres[mres]
    
    # Confidence intervals
    res$main$modCI <- bootstrapCI(fdres, B_outer = B,
      type_inner = "kfold", B_inner = 10, levels = c(.025, .975))
    
    #----- Sensitivity: model with summer temperature
    
    # Update formula
    temp_form <- paste(main_form, sprintf("Tmean_%s", region[i]), sep = " + ")
    
    # Fit model with temperature
    tempres <- FDboost(as.formula(temp_form), timeformula = NULL, 
      data = datatab, family = Gaussian(), 
      control = boost_control(nu = nu, mstop = maxstop))
    
    # Select number of boosting steps
    set.seed(1234 * i + j)
    cvres <- cvrisk(tempres, grid = 10:500,
      folds = cv(model.weights(tempres), type = "kfold", B = 10))
    mres <- mstop(cvres)
    res$temp$cv <- cvres
    res$temp$mstop <- mres
    
    # Selection of best model and kept variables
    res$temp$final_mod <- tempres[mres]
    
    #----- Output
    res
  }

# Stop parallel
stopCluster(cl)

# Save results
save.image("Results/functionalreg_results.RData")