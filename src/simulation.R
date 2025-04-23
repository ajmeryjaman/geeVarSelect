## Function for generating data

library(mvtnorm)
expit <- function(x) exp(x)/(1+exp(x))

## data.gen is a function that generates a longitudinal data set for a specific correlation
## structure. Available structures in this function are: independence, exchangeable and ar1.

# Arguments(data.gen):
#     n = Number of subjects
#     ni = A vector containing number of time points for each subject
#     sigma2.e = Error variance
#     alpha = Correlation parameter
#     corstr = The correlation structure among the repeated outcomes
#     autocorr.coef = The autocorrelation coefficient for inducing correlation among the
#     continuous confounders and the noise covariates

data.gen <- function(n, ni, sigma2.e, alpha, corstr){
  # V is the covariance matrix of the covariates (Z1, ..., Z100)
  ncovs <- 100
  Corr <- matrix(0.3, ncovs, ncovs)
  diag(Corr) <- 1
  V <- diag(rep(sqrt(0.25), ncovs))%*%Corr%*%diag(rep(sqrt(0.25), ncovs))
  beta <- c(rep(0.5, 5), rep(0.2, 5), rep(0.05, 10), rep(0, 80)) #model parameters
  
  ##outcome dependence
  corr.mat <- switch(corstr,
                     "exchangeable" = toeplitz(c(1, rep(alpha, ni-1))),
                     "ar1" = toeplitz(alpha^(0:(ni-1))),
                     "independence" = diag(ni)
  )
  cov.mat <- diag(sqrt(sigma2.e), ni) %*% corr.mat %*% diag(sqrt(sigma2.e), ni)
  
  Z <- Y <- vector(mode="list", length=n)
  for(i in 1:n){
    Z[[i]] <- rmvnorm(ni, mean = rep(0, ncovs), sigma = V)
    Y[[i]] <- rmvnorm(1, mean = Z[[i]]%*%beta, sigma = cov.mat)
  }
  Z <- do.call(rbind, Z)
  colnames(Z) <- paste("Z", 1:ncovs, sep="")
  data <- data.frame(id=rep(1:n, each=ni), Z, Y=unlist(Y))
  return(data)
}

pmse <- function(Y, X, beta){
  sum((Y - cbind(1, X)%*%beta)^2)/length(Y)
}

## Function for a single simulation and estimation for the generated data

simFUNnew <- function(n, n.i, sigma2.e, alpha.true, maxitr.penGEE, maxitr.boosting,
                      lambda.seq, working.strs, seed, ntest){
  set.seed(seed)
  data <- data.gen(n = n, ni = n.i, sigma2.e = sigma2.e, alpha = alpha.true, corstr = "exchangeable")
  
  n.wstrs <- length(working.strs)
  # mean.model <- as.formula(paste("Y ~", paste(paste("Z", 1:100, sep=""), collapse = "+"),
  #                                collapse=""))
  # Solve penalized GEE for all corr strs and under all lambda values
  out.penalizedGEE <- lapply(1:n.wstrs, function(w) lapply(lambda.seq, function(k){
    output <- penalizedGEE(Y=as.matrix(data$Y), X=as.matrix(data[,2:101]), id=data$id,
                           wc.str = working.strs[w], lambda = k, maxitr = maxitr.penGEE,
                           alpha_fixed = alpha.true)
    # print(c(w,k))
    return(output)
  }))
  
  out.GEEboosting <- lapply(1:n.wstrs, function(w){
    output <- geeBOOST(Y=as.matrix(data$Y), X=as.matrix(data[,2:101]), id=data$id,
                       wc.str = working.strs[w], alpha_fixed = alpha.true,
                       beta.init=rep(0, dim(as.matrix(data[,2:101]))[2]+1), maxitr = maxitr.boosting, eps=0.01, thresh=1)
    # print(w)
    return(output)
  })
  
  errors <- lapply(1:n.wstrs, function(w) unlist(lapply(out.penalizedGEE[[w]], function(x) x$error)))
  all.errors <- unlist(lapply(1:n.wstrs, function(w) sum(errors[[w]]==0) == 0))
  
  if(sum(all.errors)==0){
    data.test <- vector(ntest, mode="list")
    for (test in 1:ntest){
      data.test[[test]] <- data.gen(n = n, ni = n.i, sigma2.e = sigma2.e, alpha = alpha.true, corstr = "exchangeable")
    }
    
    PMSE.EEboost <- lapply(1:n.wstrs, function(w) apply(out.GEEboosting[[w]], 1, function(x)
      mean(sapply(1:ntest, function(test) 
        pmse(Y=as.matrix(data.test[[test]]$Y), X=as.matrix(data.test[[test]][,2:101]),
             beta=x)))))
    
    path.penalizedGEE <- lapply(1:n.wstrs, function(w) lapply(out.penalizedGEE[[w]], function(x){
      if(x$error==0){
        est <- x$estimate
      } else est <- NA
      return(est)
    }))
    
    PMSE.penalizedGEE <- lapply(1:n.wstrs, function(w) lapply(out.penalizedGEE[[w]], function(x){
      if(x$error==0){
        PMSE <- mean(sapply(1:ntest, function(test) 
          pmse(Y=as.matrix(data.test[[test]]$Y), X=as.matrix(data.test[[test]][,2:101]), beta=x$estimate)))
      } else PMSE <- NA
      return(PMSE)
    }))
    
    for(w in 1:n.wstrs){
      path.penalizedGEE[[w]] <- t(matrix(unlist(path.penalizedGEE[[w]]), 101, length(lambda.seq)))
      PMSE.penalizedGEE[[w]] <- unlist(PMSE.penalizedGEE[[w]])
    }
    
    all.output <- list(path.EEboost = out.GEEboosting, pmse.EEboost = PMSE.EEboost,
                       path.penalizedGEE = path.penalizedGEE, pmse.penalizedGEE = PMSE.penalizedGEE, errors = errors)
    any.error <- 0
    
  } else {
    all.output<- NULL
    any.error <- 1
  }
  return(list(all.output=all.output, any.error=any.error))
}



# simFUN <- function(n, n.i, sigma2.e, alpha.true, maxitr.penGEE, maxitr.boosting,
#                    lambda.seq, working.strs, seed){
#   set.seed(seed)
#   data <- data.gen(n = n, ni = n.i, sigma2.e = sigma2.e, alpha = alpha.true, corstr = "exchangeable")
#   
#   n.wstrs <- length(working.strs)
#   # mean.model <- as.formula(paste("Y ~", paste(paste("Z", 1:100, sep=""), collapse = "+"),
#   #                                collapse=""))
#   # Solve penalized GEE for all corr strs and under all lambda values
#   out.penalizedGEE <- lapply(1:n.wstrs, function(w) lapply(lambda.seq, function(k){
#     output <- penalizedGEE(Y=as.matrix(data$Y), X=as.matrix(data[,2:101]), id=data$id,
#                            wc.str = working.strs[w], lambda = k, maxitr = maxitr.penGEE,
#                            alpha_fixed = alpha.true)
#     print(c(w,k))
#     return(output)
#   }))
#   
#   out.GEEboosting <- lapply(1:n.wstrs, function(w){
#     output <- geeBOOST(Y=as.matrix(data$Y), X=as.matrix(data[,2:101]), id=data$id,
#                        wc.str = working.strs[w], alpha_fixed = alpha.true,
#                        beta.init=rep(0, dim(as.matrix(data[,2:101]))[2]+1), maxitr = maxitr.boosting, eps=0.01, thresh=1)
#     print(w)
#     return(output)
#   })
#   
#   errors <- lapply(1:n.wstrs, function(w) unlist(lapply(out.penalizedGEE[[w]], function(x) x$error)))
#   all.errors <- unlist(lapply(1:n.wstrs, function(w) sum(errors[[w]]==0) == 0))
#   if(sum(all.errors)==0){
#     # get the fit corresponding to the optimal lambda chosen by the information Criterion (BIC or QIC)
#     calcBIC <- lapply(1:n.wstrs, function(w) lapply(which(errors[[w]]==0), function(k)
#       BIC(out.penalized.GEE=out.penalizedGEE[[w]][[k]], Y=as.matrix(data$Y), X=as.matrix(data[,2:101]),
#           id=data$id, wc.str=working.strs[w])))
#     calcQIC1 <- lapply(1:n.wstrs, function(w) lapply(which(errors[[w]]==0), function(k)
#       QIC1(out.penalized.GEE=out.penalizedGEE[[w]][[k]], Y=as.matrix(data$Y), X=as.matrix(data[,2:101]),
#           id=data$id, wc.str=working.strs[w])))
#     calcQIC2 <- lapply(1:n.wstrs, function(w) lapply(which(errors[[w]]==0), function(k)
#       QIC2(out.penalized.GEE=out.penalizedGEE[[w]][[k]], Y=as.matrix(data$Y), X=as.matrix(data[,2:101]),
#            id=data$id, wc.str=working.strs[w])))
#     calcQIC3 <- lapply(1:n.wstrs, function(w) lapply(which(errors[[w]]==0), function(k)
#       QIC3(out.penalized.GEE=out.penalizedGEE[[w]][[k]], Y=as.matrix(data$Y), X=as.matrix(data[,2:101]),
#            id=data$id, wc.str=working.strs[w])))
#     
#     lambda.no.error <- lapply(1:n.wstrs, function(w) lambda.seq[which(errors[[w]]==0)])
#     
#     lambda.selected.QIC1 <- lapply(1:n.wstrs, function(w) lambda.no.error[[w]][which.min(unlist(calcQIC1[[w]]))])
#     lambda.selected.QIC2 <- lapply(1:n.wstrs, function(w) lambda.no.error[[w]][which.min(unlist(calcQIC2[[w]]))])
#     lambda.selected.QIC3 <- lapply(1:n.wstrs, function(w) lambda.no.error[[w]][which.min(unlist(calcQIC3[[w]]))])
#     lambda.selected.BIC <- lapply(1:n.wstrs, function(w) lambda.no.error[[w]][which.min(unlist(calcBIC[[w]]))])
# 
#     res.QIC1 <- lapply(1:n.wstrs, function(w) 
#       out.penalizedGEE[[w]][[which(lambda.seq==lambda.selected.QIC1[[w]])]])
#     res.QIC2 <- lapply(1:n.wstrs, function(w) 
#       out.penalizedGEE[[w]][[which(lambda.seq==lambda.selected.QIC2[[w]])]])
#     res.QIC3 <- lapply(1:n.wstrs, function(w) 
#       out.penalizedGEE[[w]][[which(lambda.seq==lambda.selected.QIC3[[w]])]])
#     res.BIC <- lapply(1:n.wstrs, function(w) 
#       out.penalizedGEE[[w]][[which(lambda.seq==lambda.selected.BIC[[w]])]])
#     all.output<- list(res.penalizedGEE.QIC1 = res.QIC1, res.penalizedGEE.QIC2 = res.QIC2,
#                       res.penalizedGEE.QIC3 = res.QIC3, res.penalizedGEE.BIC = res.BIC,
#                       res.GEEboosting = out.GEEboosting)
#     any.error <- 0
#   } else {
#     all.output<- NULL
#     any.error <- 1
#   }
#   return(list(all.output=all.output, any.error=any.error))
# }


