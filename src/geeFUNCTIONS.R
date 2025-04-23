#' Functions required to solve penalized GEE
#' @description This function solves penalized GEE problem for a given longitudinal data,
#'  a specific working correlation structure and a single value of the tuning parameter.
#' @param Y A vector containing response values.
#' @param X A matrix of covariates (all numeric).
#' @param id A vector containing subject identifiers.
#' @param wc.str A character string specifying the working correlation structure. Currently allowed
#'  structures: "independence", "exchangeable", "exchangeable_fixed", "ar1", "ar1_fixed", and
#'  "unstructured".
#' @param lambda A single value of the tuning parameter.
#' @param maxitr Maximum number of iterations.
#' @param penalty The penalty type to be used, available options include "SCAD" and "MCP".
#'  "SCAD" refers to the Smoothly Clipped Absolute Deviation penalty and
#'  "MCP" refers to the Minimax Concave Penalty.
#' @param alpha_fixed The known value of alpha which must be given when "exchangeable_fixed" or
#'  "ar1_fixed" is used as wc.str.
#'
#' @return A list containing the estimates, asymptotic variance etc.
#' @export


penGEE <- function(Y, X, id, wc.str, lambda, maxitr, penalty, alpha_fixed){
  # n: number of patients
  # p: number of parameters
  # ni, l, y, e: all are list of length n
  # estimate.current: a vector of length 2*p

  n <- length(unique(id))
  ni <- lapply(1:n, function(i) sum(id==id[i]))
  y <- split(Y, id)
  x <- as.data.frame(X)
  names(x) <- colnames(X)
  l.mat.split <- split(x, id)
  l <- lapply(1:n, function(i) as.matrix(cbind(1,l.mat.split[[i]])))
  p <- dim(l[[1]])[2]

  sum1 <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%y[[i]])) # p*1
  sum2 <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%l[[i]])) # p*p

  ## initial estimate, iteration = 0
  estimate.current <- solve(sum2)%*%sum1
  e <- lapply(1:n, function(i) y[[i]]-l[[i]]%*%estimate.current)
  sigma2.hat <- phi.hat.fun(n, ni, e)
  alpha.hat <- switch(wc.str,
                      "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                      "exchangeable_fixed" = alpha_fixed,
                      "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                      "ar1_fixed" = alpha_fixed,
                      "independence" = 0,
                      "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))
  itr <- 0

  ### iteration = 1

  ## corr.mat, V, V.inv >> list of length n
  V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
  V.inv <- lapply(1:n, function(i) solve(V[[i]]))
  sum1.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%e[[i]]))
  sum2.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
  penalty.der <- switch(penalty,
                        "SCAD"=q.scad(lambda, abs(estimate.current[2:p])),
                        "MCP"=q.mcp(lambda, abs(estimate.current[2:p])))
  En.mat <- diag(c(0, penalty.der/(10^(-6)+abs(estimate.current[2:p]))))
  estimate.new <- estimate.current + solve(sum2.n+n*En.mat)%*%(sum1.n-n*En.mat%*%estimate.current)
  e <- lapply(1:n, function(i) y[[i]]-l[[i]]%*%estimate.new)
  sigma2.hat <- phi.hat.fun(n, ni, e)
  alpha.hat <- switch(wc.str,
                      "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                      "exchangeable_fixed" = alpha_fixed,
                      "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                      "ar1_fixed" = alpha_fixed,
                      "independence" = 0,
                      "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))

  itr <- 0+1

  ### iterations: 2, 3, ... , till convergence/maxitr
  while(max(abs(estimate.new-estimate.current)) > 0.001 & itr<= maxitr){
    estimate.current <- estimate.new
    ## corr.mat, V, V.inv >> list of length n
    V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
    V.inv <- lapply(1:n, function(i) solve(V[[i]]))
    sum1.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%e[[i]]))
    sum2.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
    penalty.der <- switch(penalty,
                          "SCAD"=q.scad(lambda, abs(estimate.current[2:p])),
                          "MCP"=q.mcp(lambda, abs(estimate.current[2:p])))
    En.mat <- diag(c(0, penalty.der/(10^(-6)+abs(estimate.current[2:p]))))
    estimate.new <- estimate.current + solve(sum2.n+n*En.mat)%*%(sum1.n-n*En.mat%*%estimate.current)
    e <- lapply(1:n, function(i) y[[i]]-l[[i]]%*%estimate.new)
    sigma2.hat <- phi.hat.fun(n, ni, e)
    alpha.hat <- switch(wc.str,
                        "exchangeable" = corr.exch.fun(n, ni, e, phi = sigma2.hat),
                        "exchangeable_fixed" = alpha_fixed,
                        "ar1" = corr.ar1.fun(n, ni, e, phi = sigma2.hat),
                        "ar1_fixed" = alpha_fixed,
                        "independence" = 0,
                        "unstructured" = corr.un.fun(n, ni, e, phi = sigma2.hat))
    itr <- itr+1
  }

  ##Asymptotic variance calculation
  if (itr>0 & itr <= maxitr){
    V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
    V.inv <- lapply(1:n, function(i) solve(V[[i]]))
    sum2.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
    penalty.der <- switch(penalty,
                          "SCAD"=q.scad(lambda, abs(estimate.new[2:p])),
                          "MCP"=q.mcp(lambda, abs(estimate.new[2:p])))
    En.mat <- diag(c(0, penalty.der/(10^(-6)+abs(estimate.new[2:p]))))

    I.hat <- Reduce("+", lapply(1:n, function(i)
      t(l[[i]])%*%V.inv[[i]]%*%e[[i]]%*%t(e[[i]])%*%V.inv[[i]]%*%l[[i]]))
    B.hat.inv <- solve(sum2.n+n*En.mat[1:p,1:p])
    var.sand <- B.hat.inv%*%I.hat%*%B.hat.inv
    return(list(estimate = estimate.new, phi = sigma2.hat, alpha = alpha.hat,
                asymp.var = var.sand, nitr=itr, error = 0, En.mat = En.mat))
  } else {
    #print("Estimation did not converge")
    return(list(error = 1, nitr=itr))
  }
}

### Other functions to be called

## Method of moments estimate for the variance parameter sigma^2
#' Function to perform method of moments estimation for the variance parameter
#'
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param e A list of length n containing the residuals.
#'
#' @return Method of moments estimate for the variance parameter
#' @export

phi.hat.fun <- function(n, ni, e){
  sum.phi.i <- lapply(1:n, function(i) (1/ni[[i]])*sum((e[[i]])^2))
  sigma2.hat <- sum(unlist(sum.phi.i))/n
  return(sigma2.hat)
}

## Method of moments estimate for the correlation parameter alpha under exchangeable structure
#' Function to perform method of moments estimation for the correlation parameter alpha under exchangeable structure
#'
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param e A list of length n containing the residuals.
#' @param phi The estimate of variance parameter.
#'
#' @return Method of moments estimate for the correlation parameter alpha under exchangeable structure
#' @export

corr.exch.fun <- function(n, ni, e, phi){
  sum.corr.i <- lapply(1:n, function(i) (1/(ni[[i]]*(ni[[i]]-1)))*(sum(e[[i]]%*%t(e[[i]]))
                                                                   -sum(diag(e[[i]]%*%t(e[[i]])))))
  corr <- sum(unlist(sum.corr.i))/(phi*n)
  return(corr) ####
}

## Method of moments estimate for the correlation parameter alpha under AR(1) structure
#' Function to perform method of moments estimation for the correlation parameter alpha under AR(1) structure
#'
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param e A list of length n containing the residuals.
#' @param phi The estimate of variance parameter.
#'
#' @return Method of moments estimate for the correlation parameter alpha under AR(1) structure
#' @export
corr.ar1.fun <- function(n, ni, e, phi){
  sum.corr.i <- lapply(1:n, function(i) (1/(ni[[i]]-1))*(t(e[[i]][1:(ni[[i]]-1)])%*%(e[[i]][2:ni[[i]]])))
  corr <- sum(unlist(sum.corr.i))/(phi*n)
  return(corr) ####
}

## Method of moments estimate for the correlation parameters under unstructured correlation structure
#' Function to perform method of moments estimation for the correlation parameter alpha under unstructured correlation structure
#'
#' @param n A scalar value indicating the number of subjects.
#' @param ni A list of length n containing the number of measurements per subject.
#' @param e A list of length n containing the residuals.
#' @param phi The estimate of variance parameter.
#'
#' @return Method of moments estimate for the correlation parameter alpha under unstructured correlation structure
#' @export
corr.un.fun <- function(n, ni, e, phi){
  corr.i <- lapply(1:n, function(i) (1/phi)*e[[i]]%*%t(e[[i]]))
  corr.i.upper <- lapply(corr.i, function(x) cbind(x[upper.tri(x)], 1:length(x[upper.tri(x)])))
  allMat <- as.data.frame(do.call("rbind", corr.i.upper))
  names(allMat) <- c("corr", "pos")
  corr.pars <- tapply(allMat$corr, allMat$pos, mean)
  return(corr.pars)
}

## Creates correlation matrix from correlation parameters under a working correlation structure
#' Function to construct the correlation matrix for a specific subject
#'
#' @param alpha The value of correlation parameter (a scalar value for "exchangeable" or "AR(1)"
#' structure, or a vector with the elements of the upper triangular part for "unstructured").
#' @param ni A scalar value showing the number of measurements for a subject.
#' @param corstr The working correlation structure.
#'
#' @return The correlation matrix of dimension ni*ni.
#' @export

corMatF <- function(alpha, ni, corstr){
  switch(corstr,
         "exchangeable" = toeplitz(c(1, rep(alpha, ni-1))),
         "exchangeable_fixed" = toeplitz(c(1, rep(alpha, ni-1))),
         "ar1" = toeplitz(alpha^(0:(ni-1))),
         "ar1_fixed" = toeplitz(alpha^(0:(ni-1))),
         "unstructured" = vec2uMat(alpha=alpha[1:(ni*(ni-1)/2)], ni),
         "independence" = diag(ni),
  )
}

## Creates correlation matrix from correlation parameters under unstructured correlation structure
#' Function to construct correlation matrix from correlation parameters for unstructured correlation structure
#'
#' @param alpha A vector with the elements of the upper triangular part for "unstructured".
#' @param ni A scalar value showing the number of measurements for a subject.
#'
#' @return An unstructured correlation matrix
#' @export

vec2uMat <- function(alpha, ni){
  x <- matrix(1, ni, ni)
  x[upper.tri(x)] <- alpha
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  return(x)
}

### Derivative function of SCAD penalty
#' Function to construct the derivative function of SCAD penalty
#'
#' @param lambda A scalar value of the tuning parameter.
#' @param par A scalar value for a single regression parameter estimate
#'
#' @return The derivative function of SCAD penalty
#' @export

q.scad <- function(lambda, par){
  b <- 3.7
  term1 <- 1
  term2 <- ((b*lambda-par)*as.numeric((b*lambda-par)>0))/((b-1)*lambda)
  out <- lambda*(term1*as.numeric(par <= lambda)+term2*as.numeric(par > lambda))
  return(out)
}

### Derivative function of MCP penalty
#' Function to construct the derivative function of MCP penalty
#'
#' @param lambda A scalar value of the tuning parameter.
#' @param par A scalar value for a single regression parameter estimate
#'
#' @return The derivative function of MCP penalty
#' @export
q.mcp <- function(lambda, par){
  b <- 3.7
  out <- ifelse(abs(par) <= b*lambda, (lambda - abs(par)/b)*sign(par), 0)
  return(out)
}


#' Calculate BIC for a single value of lambda (tuning parameter)
#' @param out.penalized.GEE A list corresponding to the output of penalized GEE for a
#' single value of the tuning parameter. This is the output of the function penGEE().
#' @param Y A vector containing response values.
#' @param X A matrix of covariates (all numeric).
#' @param id A vector containing subject identifiers.
#' @param wc.str A working correlation structure.
#'
#' @return The value of BIC under a working correlation structure
#' @export

BIC <- function(out.penalized.GEE, Y, X, id, wc.str){
  n <- length(unique(id))
  ni <- lapply(1:n, function(i) sum(id==id[i]))
  y <- split(Y, id)
  x <- as.data.frame(X)
  names(x) <- colnames(X)
  l.mat.split <- split(x, id)
  l <- lapply(1:n, function(i) as.matrix(cbind(1,l.mat.split[[i]])))
  p <- dim(l[[1]])[2]

  estimate <- out.penalized.GEE$estimate
  sigma2.hat <- out.penalized.GEE$phi
  alpha.hat <- out.penalized.GEE$alpha
  En.mat <- out.penalized.GEE$En.mat
  V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
  V.inv <- lapply(1:n, function(i) solve(V[[i]]))
  e <- lapply(1:n, function(i) y[[i]]-l[[i]]%*%estimate)

  sum2.n <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%V.inv[[i]]%*%l[[i]]))
  df.lambda <- sum(diag(solve(sum2.n + n*En.mat)%*%sum2.n))
  sum.BIC <- Reduce("+", lapply(1:n, function(i) sum(e[[i]]^2)))
  BIC <- log(sum.BIC/sum(unlist(ni))) + log(log(n))*log(2*p)*df.lambda/n ##for high-dimension
  # BIC <- log(sum.BIC/sum(unlist(ni))) + log(n)*df.lambda/n ##works for low-dimension
  return(BIC)
}

# QIC1 <- function(out.penalized.GEE, Y, X, id, wc.str){
#   estimate <- out.penalized.GEE$estimate
#   Q <- sum((Y - cbind(1, X)%*%estimate)^2/(-2))
#   qic1 <- -2*Q + 2*(sum(estimate!=0)-1)
#   return(qic1)
# }


#' Calculate QIC for a single value of lambda (tuning parameter)
#' @param out.penalized.GEE A list corresponding to the output of penalized GEE for a
#' single value of the tuning parameter. This is the output of the function penGEE().
#' @param Y A vector containing response values.
#' @param X A matrix of covariates (all numeric).
#' @param id A vector containing subject identifiers.
#' @param wc.str A working correlation structure.
#'
#' @return The value of QIC under a working correlation structure
#' @export

QIC <- function(out.penalized.GEE, Y, X, id, wc.str){
  n <- length(unique(id))
  x <- as.data.frame(X)
  names(x) <- colnames(X)
  l.mat.split <- split(x, id)
  l <- lapply(1:n, function(i) as.matrix(cbind(1,l.mat.split[[i]])))
  estimate <- out.penalized.GEE$estimate
  sigma2.hat <- out.penalized.GEE$phi
  Q <- sum((Y - cbind(1, X)%*%estimate)^2/(-2))
  sum.qic <- Reduce("+", lapply(1:n, function(i) t(l[[i]])%*%l[[i]]))/sigma2.hat
  qic.pan <- -2*Q + 2*sum(diag(sum.qic%*%out.penalized.GEE$asymp.var))
  return(qic.pan)
}
