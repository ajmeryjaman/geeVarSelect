rm(list = ls())

sensitivity <- function(beta.est){
  beta.true <- c(rep(0.5, 5), rep(0.2, 5), rep(0.05, 10), rep(0, 80))
  return(sum((beta.true != 0)*(beta.est != 0))/sum(beta.true != 0))
}

specificity <- function(beta.est){
  beta.true <- c(rep(0.5, 5), rep(0.2, 5), rep(0.05, 10), rep(0, 80))
  return(sum((beta.true == 0)*(beta.est == 0))/sum(beta.true == 0))
}

summaryFN <- function(out, alpha.true){
  nsim <- 100
  res.sensitivityIN <- list(EEboost = rep(NA, nsim), penalized = rep(NA, nsim))
  res.specificityIN <- list(EEboost = rep(NA, nsim), penalized = rep(NA, nsim))
  res.sensitivityEF <- list(EEboost = rep(NA, nsim), penalized = rep(NA, nsim))
  res.specificityEF <- list(EEboost = rep(NA, nsim), penalized = rep(NA, nsim))
  res.sensitivityEE <- list(EEboost = rep(NA, nsim), penalized = rep(NA, nsim))
  res.specificityEE <- list(EEboost = rep(NA, nsim), penalized = rep(NA, nsim))

  for(sim in 1:nsim){

    path <- out[[sim]]$all.output$path.EEboost[[1]]
    pmse.val <- out[[sim]]$all.output$pmse.EEboost[[1]]
    res.sensitivityIN$EEboost[sim] <- sensitivity(beta.est = path[which.min(pmse.val),2:101])
    res.specificityIN$EEboost[sim] <- specificity(beta.est = path[which.min(pmse.val),2:101])

    path <- out[[sim]]$all.output$path.penalizedGEE[[1]]
    pmse.val <- out[[sim]]$all.output$pmse.penalizedGEE[[1]]
    est.pen <- path[which.min(pmse.val),2:101]
    est.pen[abs(est.pen) <= 0.001] <- 0
    res.sensitivityIN$penalized[sim] <- sensitivity(beta.est = est.pen)
    res.specificityIN$penalized[sim] <- specificity(beta.est = est.pen)

    path <- out[[sim]]$all.output$path.EEboost[[2]]
    pmse.val <- out[[sim]]$all.output$pmse.EEboost[[2]]
    res.sensitivityEF$EEboost[sim] <- sensitivity(beta.est = path[which.min(pmse.val),2:101])
    res.specificityEF$EEboost[sim] <- specificity(beta.est = path[which.min(pmse.val),2:101])

    path <- out[[sim]]$all.output$path.penalizedGEE[[2]]
    pmse.val <- out[[sim]]$all.output$pmse.penalizedGEE[[2]]
    est.pen <- path[which.min(pmse.val),2:101]
    est.pen[abs(est.pen) <= 0.001] <- 0
    res.sensitivityEF$penalized[sim] <- sensitivity(beta.est = est.pen)
    res.specificityEF$penalized[sim] <- specificity(beta.est = est.pen)

    path <- out[[sim]]$all.output$path.EEboost[[3]]
    pmse.val <- out[[sim]]$all.output$pmse.EEboost[[3]]
    res.sensitivityEE$EEboost[sim] <- sensitivity(beta.est = path[which.min(pmse.val),2:101])
    res.specificityEE$EEboost[sim] <- specificity(beta.est = path[which.min(pmse.val),2:101])

    path <- out[[sim]]$all.output$path.penalizedGEE[[3]]
    pmse.val <- out[[sim]]$all.output$pmse.penalizedGEE[[3]]
    est.pen <- path[which.min(pmse.val),2:101]
    est.pen[abs(est.pen) <= 0.001] <- 0
    res.sensitivityEE$penalized[sim] <- sensitivity(beta.est = est.pen)
    res.specificityEE$penalized[sim] <- specificity(beta.est = est.pen)

    print(sim)
  }
  return(list(res.sensitivityIN=res.sensitivityIN, res.specificityIN=res.specificityIN,
              res.sensitivityEF=res.sensitivityEF, res.specificityEF=res.specificityEF,
              res.sensitivityEE=res.sensitivityEE, res.specificityEE=res.specificityEE))
}

rm(out)
load("results/out_rho0.Rdata")
apply(out[[100]]$all.output$path.EEboost[[1]], 1, function(x) sum(abs(x)))
res0 <- summaryFN(out = out, alpha.true = 0)

rm(out)
load("results/out_rho03.Rdata")
res03 <- summaryFN(out = out, alpha.true = 0.3)

rm(out)
load("results/out_rho07.Rdata")
res07 <- summaryFN(out = out, alpha.true = 0.7)

rm(out)
load("results/out_rho09.Rdata")
res09 <- summaryFN(out = out, alpha.true = 0.9)

#######################
####### Table 2 ######
#######################

sens_spec <- unlist(lapply(res0, function(x) lapply(x, function(xx) median(xx))))
sens <- sens_spec[c(1,2,5,6,9,10)]
spec <- sens_spec[c(3,4,7,8,11,12)]

sens_spec <- unlist(lapply(res03, function(x) lapply(x, function(xx) median(xx))))
sens <- rbind(sens, sens_spec[c(1,2,5,6,9,10)])
spec <- rbind(spec, sens_spec[c(3,4,7,8,11,12)])

sens_spec <- unlist(lapply(res07, function(x) lapply(x, function(xx) median(xx))))
sens <- rbind(sens, sens_spec[c(1,2,5,6,9,10)])
spec <- rbind(spec, sens_spec[c(3,4,7,8,11,12)])

sens_spec <- unlist(lapply(res09, function(x) lapply(x, function(xx) median(xx))))
sens <- rbind(sens, sens_spec[c(1,2,5,6,9,10)])
spec <- rbind(spec, sens_spec[c(3,4,7,8,11,12)])


colnames(sens) <- c("IN.EE", "IN.PE", "EF.EE", "EF.PE", "EE.EE", "EE.PE")
colnames(spec) <- c("IN.EE", "IN.PE", "EF.EE", "EF.PE", "EE.EE", "EE.PE")

# my results
round(sens, 2)
round(spec, 2)

# original results
sens.p <- cbind(rep(0.55, 4), c(0.85, 0.8, 0.75, 0.8), c(0.55, 0.55, 0.65, 0.75),
                c(0.85, 0.8, 0.8, 0.85), c(0.55, 0.55, 0.60, 0.75), c(0.85, 0.8, 0.8, 0.8))
spec.p <- cbind(c(0.8, 0.83, 0.83, 0.83), c(0.37, 0.4, 0.4, 0.4), c(0.8, 0.83, 0.77, 0.77),
                c(0.37, 0.43, 0.43, 0.4), c(0.8, 0.83, 0.8, 0.77), c(0.37, 0.4, 0.4, 0.4))


library(xtable)

print(xtable(cbind(c(0,0.3,0.7,0.9), sens.p[,1], sens[,1], NA, sens.p[,2], sens[,2], NA,
                   sens.p[,3], sens[,3], NA,sens.p[,4], sens[,4], NA, sens.p[,5], sens[,5], NA,
                   sens.p[,6], sens[,6]), digits = c(0,1,rep(2,17))), include.rownames = FALSE)


print(xtable(cbind(c(0,0.3,0.7,0.9), spec.p[,1], spec[,1], NA, spec.p[,2], spec[,2], NA, spec.p[,3], spec[,3], NA,
                   spec.p[,4], spec[,4], NA, spec.p[,5], spec[,5], NA, spec.p[,6], spec[,6]), digits = c(0,1,rep(2,17))), include.rownames = FALSE)








