#' UBayFS feature selection
#' @description genetic algorithm to train UBayFS feature selection model
#' @param x a UBaymodel created by build.UBaymodel
#' @return a UBaymodel with an additional list element output containing the maximum a-posteriori estimate (map)
#' @importFrom GA ga
#' @importFrom DirichletReg ddirichlet
#' @export train

train <- function(x){
  UseMethod("train")
}

#' @method train UBaymodel
#' @export

train.UBaymodel = function(x){

  if(class(x) != "UBaymodel"){
    stop("Wrong class of model")
  }

  # define prior weights
  alpha = as.numeric(x$user.params$weights)

  # define ensemble counts
  #delta = as.numeric(x$ensemble.params$output)

  if(x$optim.params$method == "GA"){
    print("Running Genetic Algorithm")
    x$output <- train_GA(alpha,
                         #delta,
                         x$ensemble.params$output,
                         x$user.params$constraints,
                         x$user.params$block_constraints,
                         x$optim.params,
                         colnames(x$data))
  }
  else if(x$optim.params$method == "MH"){
    print("Running Metropolis-Hastings Algorithm")
    x$output <- train_MH(alpha,
                         #delta,
                         x$ensemble.params$output,
                         x$user.params$constraints,
                         x$user.params$block_constraints,
                         x$optim.params,
                         colnames(x$data))
  }
  else{
    stop("Error: method not supported.")
  }

  return(x)
}

train_GA <- function(weights, counts, constraints, block_constraints, optim_params, feat_names){

  target_fct <- function(x){
    if(length(x) == length(weights)){
      state = x
      covariates = NULL
    }
    else{
      state = x[1:length(weights)]
      covariates = matrix(NA, nrow = length(weights), ncol = length(weights))
      covariates[upper.tri(covariates)] <- x[(length(weights) + 1): length(x)]
    }
    return(posterior(state,
                     constraints = constraints,
                     block_constraints = block_constraints,
                     weights = weights,
                     counts = counts,
                     covariates = covariates,
                     log = TRUE))}

  # Greedy algorithm to select starting vectors
  x_start = sampleInitial(post_scores = weights + colSums(counts),
                          constraints = constraints,
                          block_constraints = block_constraints,
                          size = optim_params$popsize)

  optim = ga(type = "binary",								# use GA for optimization
             fitness = target_fct,
             lower = 0,
             upper = 1,
             nBits = length(weights) + length(weights) * (length(weights) - 1)/2,
             maxiter = optim_params$maxiter,
             popSize = optim_params$popsize,
             suggestions = x_start
  )
  x_optim = optim@solution									# extract solution
  if(is.vector(x_optim)){
    x_optim <- as.data.frame(t(x_optim))
  }
  else{
    x_optim <- as.data.frame(x_optim)
  }
  colnames(x_optim) <- feat_names

  return(list(map = x_optim))
}

#' @importFrom plyr count

train_MH = function(weights, counts, constraints, block_constraints, optim_params, feat_names){

  target_fct <- function(x){
    if(length(x) == length(weights)){
      state = x
      covariates = NULL
    }
    else{
      state = x[1:length(weights)]
      covariates = matrix(NA, nrow = length(weights), ncol = length(weights))
      covariates[upper.tri(covariates)] <- x[(length(weights) + 1): length(x)]
    }
    return(posterior(state,
                     constraints = constraints,
                     block_constraints = block_constraints,
                     weights = weights,
                     counts = counts,
                     covariates = covariates,
                     log = TRUE))}

  # empirical density of proposal
  proposal_sample_size <- optim_params$popsize * optim_params$maxiter
  proposal_sample <- sampleInitial(post_scores = weights + colSums(counts),
                                    constraints = constraints,
                                    block_constraints = block_constraints,
                                    size = proposal_sample_size)
  proposal_sample <- as.data.frame(proposal_sample)
  proposal_sample$str <- apply(proposal_sample, 1, paste, collapse = ",")
  proposal_sample <- plyr::count(proposal_sample)
  g <- function(state, eps = 0.1 / proposal_sample_size, log = TRUE){
    state_str <- paste(state, collapse = ",")
    i <- which(proposal_sample$str == state_str)
    if(length(i) > 0){
      return(log(proposal_sample$freq[i] / proposal_sample_size + eps))
    }
    else{
      return(log(eps))
    }
  }


  # sample from proposal density
  X <- matrix(, nrow = 0, ncol = length(weights))        # MCMC history
  x_t = sampleInitial(post_scores = weights + colSums(counts),
                      constraints = constraints,
                      block_constraints = block_constraints,
                      size = optim_params$popsize)
  f_x_t = apply(x_t, 1, target_fct)
  g_x_t = apply(x_t, 1, g)
  acceptances <- 0

  # calculate MH-ratio
  for(t in 1:optim_params$maxiter){
    print(t)
    # new sample from proposal density
    x_new = sampleInitial(post_scores = weights + colSums(counts),
                          constraints = constraints,
                          block_constraints = block_constraints,
                          size = optim_params$popsize)
    # calculate MH-ratio
    f_x_new = apply(x_new, 1, target_fct)
    g_x_new = apply(x_new, 1, g)
    mh_ratio = apply(cbind(f_x_new - f_x_t + g_x_t - g_x_new, 0), 1, min)
    mh_ratio = exp(mh_ratio)
    acceptance = apply(cbind(mh_ratio, 1-mh_ratio), 1, sample, x = c(1,0), size = 1, replace = FALSE) == 1
    acceptances <- c(acceptances, mean(acceptance))
    print(sum(acceptance))
    x_t[acceptance,] <- x_new[acceptance,]
    X <- rbind(X, x_t)
    f_x_t[acceptance] <- f_x_new[acceptance]
  }

  # prepare output
  colnames(X) <- feat_names
  unique_X <- as.data.frame(X)
  unique_X <- plyr::count(unique_X)
  max_freq <- max(unique_X$freq)
  unique_X <- unique_X[unique_X$freq == max_freq, -ncol(unique_X)]
  colnames(unique_X) <- c(feat_names)

  return(list(post.sample = X,
              map = unique_X,
              acceptance_percent = mean(acceptances)))

}
