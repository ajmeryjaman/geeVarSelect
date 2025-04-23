#' Main function to solve penalized GEE 
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

## The main function solving penalized GEE for a given data, working
## correlation structure and a sequence of tuning parameter values

penalizedGEE <- function(Y, X, id, wc.str, lambda, maxitr = 100, penalty = "SCAD",
                         alpha_fixed = NULL){
  if(!wc.str %in% c("independence", "exchangeable", "ar1", "unstructured",
                    "exchangeable_fixed", "ar1_fixed")){
    stop("Working correlation structure must be one among independence,
         exchangeable, ar1, and unstructured\n")
  }
  data <- data.frame(Y, X, id)
  if(any(c(is.na(data)))){
    stop("Data can not contain any missing values\n")
  }
  
  # X <- apply(X, 2, function(i) (i-mean(i))/sd(i))
  
  ##Next we solve penalized GEE 
  out.penGEE <- penGEE(Y=Y, X=X, id=id, wc.str = wc.str, lambda = lambda, maxitr = maxitr, penalty = penalty,
                  alpha_fixed = alpha_fixed)
  return(out.penGEE)
}

##' Penalized GEE
##' @description This function solves penalized GEE problem for a given longitudinal data,
##' a specific working correlation structure and a sequence of tuning parameters, and returns result
##' under the optimal value (selected by a BIC-type  criterion) of the tuning parameter in the given range.
##' Currently, the function allows a continuous outcome.
##' @param Y A vector containing response values.
##' @param X A matrix of covariates (all numeric).
##' @param id A vector containing subject identifiers.
##' @param wc.str A character string specifying the working correlation structure. The
##' following are currently allowed: "independence", "exchangeable", "ar1", and "unstructured".
##' @param lambda.seq A sequence of values (in decreasing order) for the tuning parameter. Typically,
##' the first element in this sequence is a lowest positive number such that all the effect modifiers (EMs) are eliminated by the method, while
##' the last element is a value close to zero for which none of the EMs are eliminated.
##' @param maxitr Maximum number of iterations allowed. The default value is 100.
##' @param penalty The penalty type to be used, available options include "SCAD" and "MCP". The default choice is "SCAD".
##' "SCAD" refers to the Smoothly Clipped Absolute Deviation penalty and "MCP" refers to the Minimax Concave Penalty.
##' 
##' @return A list containing the following:
##' \item{estimate}{The vector of parameter estimates. First half corresponds to the blip coefficients
##' and the second half corresponds to the coefficients of the treatment-free model}
##' \item{Selected.covs}{A vector showing which variables are selected}
##' \item{sigma2.hat}{The estimated variance parameter sigma^2.}
##' \item{alpha.hat}{The estimated correlation parameter(s) alpha(s) if the provided structure is either
##' "exchangeable", "ar1, or "unstructured". For unstructured, the elements of alpha.hat correspond
##' to the upper triangular portion of working correlation matrix having dimension equal to
##' the largest cluster size.}
##' \item{asymp.var}{The sandwich variance-covariance matrix of the blip coefficients. Although, the
##' function will provided
##' sandwich variance, a data analyst should note that tests or confidence intervals based on this
##' sandwich estimates are expected to exhibit
##' inflated type I errors in finite samples.}
##' \item{nitr}{The number of iterations at which the estimation converged with the optimal tuning parameter.}
##' \item{lambda.optimal}{The optimal tuning parameter from the given sequence (lambda.seq).}
##' 
##' @export
##' 
##' @importFrom stats model.matrix as.formula binomial glm predict.glm toeplitz
##' @examples


# penalizedGEE <- function(Y, X, id, wc.str, lambda.seq, maxitr = 100, penalty = "SCAD"){
#   if(!wc.str %in% c("independence", "exchangeable", "ar1", "unstructured")){
#     stop("Working correlation structure must be one among independence,
#          exchangeable, ar1, and unstructured\n")
#   }
#   data <- data.frame(Y, X, id)
#   if(any(c(is.na(data)))){
#     stop("Data can not contain any missing values\n")
#   }
# 
#   ##Next we solve penalized GEE for a sequence of tuning parameters and
#   ##we record if there is any error (i.e., the estimation did not converge)
#   out.penGEE <- lapply(lambda.seq, function(k){
#     pG <- penGEE(Y=Y, X=X, id=id, wc.str = wc.str, lambda = k, maxitr = maxitr, penalty = penalty)
#     print(k)
#     return(pG)
#     })
#   errors <- unlist(lapply(out.penGEE, function(x) x$error))
# 
#   ##Next we split the data and construct required quantities as a list of length n
#   ## for computing the values of BIC
#   # n <- length(unique(id))
#   # ni <- lapply(1:n, function(i) sum(id==id[i]))
#   # y <- split(Y, id)
#   # x <- as.data.frame(X)
#   # names(x) <- colnames(X)
#   # l.mat.split <- split(x, id)
#   # l <- lapply(1:n, function(i) as.matrix(cbind(1,l.mat.split[[i]])))
#   # p <- dim(l[[1]])[2]
# 
#   ###BIC calculation for the tuning parameters for which the estimation converged (i.e., error = 0)
#   # BIC <- lapply(which(errors == 0), function(k)
#   #   calcBIC(out.penGEE[[k]], n, ni, p, y, l, wc.str = wc.str)
#   # )
#   
#   QIC <- lapply(which(errors == 0), function(k)
#     QIC.fn(Y, cbind(1, X), out.penGEE[[k]]$estimate)
#   )
#   lambda.no.error <- lambda.seq[which(errors == 0)]
#   lambda.selected <- lambda.no.error[which.min(unlist(QIC))]
#   res <- out.penGEE[[which(lambda.seq==lambda.selected)]]
# 
#   estimate <- res$estimate
#   sigma2.hat <- res$phi
#   alpha.hat <- ifelse(wc.str %in% c("exchangeable", "ar1", "unstructured"), res$alpha, NA)
#   asymp.var <- res$asymp.var
#   nitr <- res$nitr
# 
#   row.names(estimate) <- c("Intercept", colnames(X))
#   row.names(asymp.var) <- colnames(asymp.var) <- c("Intercept", colnames(X))
#   selected.covs <- colnames(X)[abs(res$estimate[2:p]) > 0.001]
# 
#   return(list(estimate = estimate, selected.covs = selected.covs, sigma2.hat = sigma2.hat,
#               alpha.hat = alpha.hat, asymp.var = asymp.var,
#               nitr=nitr, lambda.optimal = lambda.selected))
# }


