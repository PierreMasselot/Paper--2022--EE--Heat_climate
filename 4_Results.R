#################################################################
#
#           AUDACE - Summer mortality prediction
#                     Results summary
#
#################################################################

load("Results/functionalreg_results.RData")

#----------------------------
# Tables with AFs
#----------------------------

# Loop on regions
for (i in seq_along(region)){
  
  # Get AF variables
  afvars <- datatab[grep(sprintf("AF_.*_%s", region[i]), 
    names(datatab), value = T)]
  
  # Pretty number for AF with CIs
  afpretty <- sapply(seq_along(afvars), function(j){
    sprintf("%2.1f (%2.1f - %2.1f)", afvars[[j]] * 100,
      afcis[[i]][1,,j] * 100, afcis[[i]][2,,j] * 100)
  })
  
  # Add year
  afpretty <- cbind(years, afpretty)
  
  # Names
  colnames(afpretty) <- c("Year", sapply(strsplit(names(afvars), "_"), "[", 2))
  
  # Output
  write.table(afpretty, file = sprintf("Figures/allAF_%s.csv", region[i]),
    sep = ",", row.names = F, quote = F)
}

#----------------------------
# Score fitting
#----------------------------

#----- Cross-validated R2

# R2 function
r2fun <- function(object) 1 - (var(object$resid()) / var(object$response))

# Set seed
set.seed(4321)

# Compute CV R2
cvR2 <- lapply(reslist, sapply, function(x) {
  cvr2 <- cvrisk(x$main$final_mod, grid = mstop(x$main$final_mod), fun = r2fun,
    folds = cv(model.weights(x$main$final_mod), type = "kfold", B = 10))
  c(mean(unlist(cvr2)), sd(unlist(cvr2)))
})

#----- R2 of model with mean temperature

tempR2 <- lapply(reslist, sapply, function(x) {
  cvr2 <- cvrisk(x$temp$final_mod, grid = mstop(x$temp$final_mod), fun = r2fun,
    folds = cv(model.weights(x$temp$final_mod), type = "kfold", B = 10))
  c(mean(unlist(cvr2)), sd(unlist(cvr2)))
})

#----- Mstop of each model

mstops <- sapply(reslist, sapply, function(x) x$main$mstop)

#----------------------------
# Selected variables for each model
#----------------------------

#----- Main model

# Index of selected variables
selinds <- lapply(reslist, lapply, 
  function(x) x$main$final_mod$which(usedonly = T))

# Variables
sel_vars <- lapply(selinds, lapply, function(x){
  c(climate_labs, "Time")[x]
})

#----- Selected variables for model with temperature

# Index of selected variables
selinds <- lapply(reslist, lapply, 
  function(x) x$temp$final_mod$which(usedonly = T))

# Variables
sel_vars_temp <- lapply(selinds, lapply, function(x){
  c(climate_labs, "Time", "Temperature")[x]
})

#----------------------------
# Table 1: Summary statistics
#----------------------------

sum_tab <- sapply(dlist, function(d){
  c(totdeaths = sum(d$Death), summary(d$Tmean))
})

formatC(sum_tab, format = "f")

AFlist<- sapply(datatab[grep("AF", names(datatab))], function(x){
  afsum <- summary(x) * 100
  sprintf("%1.1f (%1.1f - %1.1f)", afsum["Mean"], afsum["Min."], afsum["Max."])
})

AFtab <- matrix(AFlist, ncol = 2, 
  dimnames = list(c("MMT", heatper * 100), region))

sapply(reslist, sapply, function(x) mstop(x$final_mod))
