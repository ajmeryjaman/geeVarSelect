#' Function to perform GEE boosting
#' @description The function performs GEE boosting algorithm for a given data, working correlation structure
#' @param Y A vector containing response values.
#' @param X A matrix of covariates (all numeric).
#' @param id A vector containing subject identifiers.
#' @param wc.str A character string specifying the working correlation structure. Currently allowed
#'  structures: "independence", "exchangeable", "exchangeable_fixed", "ar1", "ar1_fixed", and
#'  "unstructured".
#' @param alpha_fixed The known value of alpha which must be given when "exchangeable_fixed" or
#'  "ar1_fixed" is used as wc.str.
#' @param beta.init Initial parameter values.
#' @param maxitr Maximum number of iterations.
#' @param eps The value for epsilon-boosting.
#' @param thresh Threshold value for thresholded boosting.
#'
#' @return A list containing the estimates, asymptotic variance etc.
#' @export
#' 

geeBOOST <- function(Y, X, id, wc.str, alpha_fixed = NULL, beta.init, maxitr, eps, thresh){
  # n: number of patients
  # p: number of parameters
  # ni, l, y, e: all are list of length n
  # estimate.current: a vector of length 2*p
  
  n <- length(unique(id))
  ni <- lapply(1:n, function(i) sum(id==id[i]))
  y <- split(Y, id)
  x <- as.data.frame(apply(X, 2, function(i) (i-mean(i))/sd(i))) ## standardizing the covariate matrix
  names(x) <- colnames(X)
  l.mat.split <- split(x, id)
  l <- lapply(1:n, function(i) as.matrix(cbind(1,l.mat.split[[i]])))
  p <- dim(l[[1]])[2]
  
  estimate <- matrix(NA, maxitr, length(beta.init))
  
  ## initial estimate, iteration = 1
  estimate.current <- beta.init
  e <- lapply(1:n, function(i) y[[i]]-l[[i]]%*%estimate.current)
  sigma2.hat <- phi.hat.fun(n, ni, e)
  alpha.hat <- switch(wc.str,
                      "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                      "exchangeable_fixed" = alpha_fixed,
                      "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                      "ar1_fixed" = alpha_fixed,
                      "independence" = 0,
                      "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))
  itr <- 1
  
  ### iteration = 2
  V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
  V.inv <- lapply(1:n, function(i) solve(V[[i]]))
  g.n <- round(Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%e[[i]])), 8)
  
  ## why added randomizer runif, and why used threshold: threshold*max
  #estimate.new <- estimate.current + eps*sign(g.n)*(abs(g.n) >= max(abs(g.n)))
  estimate.new <- estimate.current + runif(1,0.5,1.5)*eps*sign(g.n)*(abs(g.n) >= thresh*max(abs(g.n)))
  e <- lapply(1:n, function(i) y[[i]]-l[[i]]%*%estimate.new)
  sigma2.hat <- phi.hat.fun(n, ni, e)
  alpha.hat <- switch(wc.str,
                      "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                      "exchangeable_fixed" = alpha_fixed,
                      "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                      "ar1_fixed" = alpha_fixed,
                      "independence" = 0,
                      "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))
  
  estimate[itr,] <- estimate.new
  itr <- itr+1
  
  ### iterations: 3, 4, ...
  while(itr <= maxitr){
    estimate.current <- estimate.new
    V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
    V.inv <- lapply(1:n, function(i) solve(V[[i]]))
    g.n <- round(Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%e[[i]])), 8)
    
    #estimate.new <- estimate.current + eps*sign(g.n)*(abs(g.n) >= max(abs(g.n)))
    estimate.new <- estimate.current + runif(1,0.5,1.5)*eps*sign(g.n)*(abs(g.n) >= thresh*max(abs(g.n)))
    e <- lapply(1:n, function(i) y[[i]]-l[[i]]%*%estimate.new)
    sigma2.hat <- phi.hat.fun(n, ni, e)
    alpha.hat <- switch(wc.str,
                        "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                        "exchangeable_fixed" = alpha_fixed,
                        "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                        "ar1_fixed" = alpha_fixed,
                        "independence" = 0,
                        "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))
    estimate[itr, ] <- estimate.new
    # print(itr)
    itr <- itr+1
  }
  return(estimate)
}




