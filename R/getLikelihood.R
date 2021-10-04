#' likelihood
#' @description evaluates the value of the model likelihood
#' @param state a binary membership vector describing a feature set
#' @param param a vector of posterior weights for the Multinomial function
#' @param log whether the admissibility should be returned on log scale
#' @return a posterior probability value
#' @importFrom stats dmultinom
#' @importFrom MM MM_single
#' @importFrom MM paras
#' @importFrom MM p
#' @importFrom MM theta

likelihood = function(state, param, log = TRUE){								# target function for optimization procedure
  if(!is.matrix(state)){
    stop("Error: state must be a matrix")
  }
  if(!is.list(param) || length(param) != 2){
    stop("Error: param must be a list with 2 elements")
  }
  n = length(param[[1]])
  if(is.null(param[[2]])){
    lik <- apply(state, 1, function(s){return(dmultinom(s, sum(s), prob = param[[1]], log = TRUE))})
    lik <- sum(lik)
  }
  else{
    par <- MM::paras(n)
    MM::p(par) <- param[[1]]/sum(param[[1]])
    MM::theta(par) <- param[[2]]
    lik <- apply(state, 1, MM_single, paras = par, givelog = TRUE)
    lik <- sum(lik)
  }
  if(!log){
    lik <- exp(lik)
  }
  return(lik)
}

#' Model likelihood
#' @description evaluates the value of the model likelihood
#' @param model a UBaymodel object created using build.UBaymodel
#' @param state a binary membership vector describing a feature set
#' @param log whether the admissibility should be returned on log scale
#' @return a posterior probability value
#' @export

getLikelihood <- function(state, covariates = NULL, model, log = TRUE){
  return(likelihood(
    state = model$ensemble.params$output,
    param = list(state, covariates),
    log = log
  ))
}
