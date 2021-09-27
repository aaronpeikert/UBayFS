#' likelihood
#' @description evaluates the value of the model likelihood
#' @param state a binary membership vector describing a feature set
#' @param param a vector of posterior weights for the Multinomial function
#' @param dmax an positive integer denoting the maximum degree of feature interactions considered in the model
#' @param log whether the admissibility should be returned on log scale
#' @return a posterior probability value
#' @importFrom stats dmultinom

likelihood = function(state, param, dmax, log = TRUE){								# target function for optimization procedure
  n = length(param)
  lik <- 0
  for(d in 1:dmax){
    n = choose(length(param), d)

    if(d > 1){
      trans_mat = matrix(0, nrow = n, ncol = length(param))

      ind_mat <- t(combn(1:length(param), d))

      trans_mat[cbind(1:n, as.vector(ind_mat))] <- 1
      theta_d = as.vector(param %*% t(trans_mat))

    }
    else{
      theta_d = as.vector(param)
    }
    lik <- lik + dmultinom(state[[d]], size = sum(state[[d]]), prob = theta_d, log = TRUE)
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

getLikelihood <- function(state, model, log = TRUE){
  return(likelihood(
    state = model$ensemble.params$output,
    param = state,
    dmax = model$ensemble.params$input$dmax,
    log = log
  ))
}
