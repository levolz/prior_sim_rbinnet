# build_IM() for Ising model graph ------------------------------------
build_IM <- function(p, prob){
  # create true model graph matrix
  Graph <- matrix(data = 0, nrow = p, ncol = p)
  Graph[lower.tri(Graph)] <- rbinom(n = p * (p - 1) / 2,
                                    size = 1,
                                    prob = prob)
  Graph <- Graph * runif(n = p ^ 2, min = 0.5, max = 2)
  Graph <- Graph + t(Graph)
  # sum over edges to get node thresholds
  Thresholds <- -rowSums(Graph) / 2
  
  return(list(Graph = Graph, Thresholds = Thresholds))
}

# Ising Estimation function -------------------------------------------
# output - estimated parameters, their standard deviation, inclusions probability
estimate_IM <- function(Data, MCMC_iter = 1e5, precision = .975){
  fit_screenAndSelect <- rbinnet::select_structure(x = Data,                            #(M) n * p matrix of binary responses
                                                   sigma = sigma.map,                   #(O) p * p matrix of Ising parameters
                                                   theta = 0.5,                         #(O) prior inclusion probability
                                                   prior_var_intercepts = 1,
                                                   #slab_var = nu1,                     #(M) p * p matrix; prior slab variance
                                                   #spike_var = nu0,                    #(M) p * p matrix; prior spike variance
                                                   #precision = c0.GM,                  #(M) penalty parameter
                                                   output_samples = T,                  #(M) if TRUE, outputs posterior draws
                                                   number_iterations = MCMC_iter,      #(M) no. iterations Gibbs sampler
                                                   number_burnin_iterations = 0,        #(M) no. burnin iterations Gibbs sampler
                                                   #include = fit_screen$gamma > 0.5,   #(O) p * p matrix: 1 (0) = include (exclude) in analysis
                                                   hierarchical = FALSE, 
                                                   precision = precision)
  
  # Parameters of the model (interaction parameters on the off-diagonal; thresholds on the diagonal)
  parameters <- vector_to_matrix(colMeans(fit_screenAndSelect$sigma_samples), 
                                 p, diag = T)
  param_sd <- vector_to_matrix(colSD(fit_screenAndSelect$sigma_samples), 
                               p, diag = T)
  
  #Inclusion probability per edges
  inc_probs_screenAndSelect <- vector_to_matrix(fit_screenAndSelect$posterior_probability %*% fit_screenAndSelect$structures, p)
  
  ## combine model output in list
  model_output <- list(parameters = parameters,
                       parameters_stdev = param_sd,
                       inclusion_prob = inc_probs_screenAndSelect)
  return(model_output)
}
