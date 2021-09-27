#' posterior function
#' @description evaluates the value of the posterior probability density (target function)
#' @param state a binary membership vector describing a feature set
#' @param constraints a list constaining a matrix A and a vector b representing the inequality system Ax<=b and a vector rho
#' @param block_constraints a list containing a relaxed system Ax<=b of user constraints on feature blocks, given as matrix A, vector b and vector or scalar rho (relaxation parameters); see buildConstraints function
#' @param weights a vector of prior weights for the Dirichlet distribution
#' @param counts a vector of likelihood counts for the Multinomial distribution
#' @param dmax an positive integer denoting the maximum degree of feature interactions considered in the model
#' @param log whether the admissibility should be returned on log scale
#' @return a posterior probability value
#' @importFrom DirichletReg ddirichlet

posterior = function(state, constraints, block_constraints, weights, counts, dmax, log = TRUE){								# target function for optimization procedure
  post <-
    admissibility(state, 									# log-admissibility function
                  constraints,
                  sum(weights + counts[[1]]),
                  log = TRUE) +
    block_admissibility(state, 									# log-admissibility function
                        block_constraints,
                        sum(weights + counts[[1]]) / nrow(block_constraints$block_matrix),
                        log = TRUE) +
    ddirichlet(t(state + 0.01), 							# log-dirichlet-density (with small epsilon to avoid errors from 0 probs)
                 alpha = weights,
                 log = TRUE)+
    likelihood(state = counts,
               param = t(state + 0.01),
               dmax = dmax,
               log = TRUE)  # likelihood
  if(!log){
    post <- exp(post)
  }
  return(post)
}

#' Model posterior
#' @description evaluates the value of the posterior probability density (target function)
#' @param model a UBaymodel object created using build.UBaymodel
#' @param state a binary membership vector describing a feature set
#' @param log whether the admissibility should be returned on log scale
#' @return a posterior probability value
#' @export

getPosterior <- function(state, model, log = TRUE){
  return(posterior(
    state = state,
    constraints = model$user.params$constraints,
    block_constraints = model$user.params$block_constraints,
    weights = model$user.params$weights,
    counts = model$ensemble.params$output,
    dmax = model$ensemble.params$input$dmax,
    log = log
  ))
}
