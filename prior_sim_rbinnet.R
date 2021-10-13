# Load Packages -------------------------------------------------------
library("BGGM")
library("BDgraph")
library("rbinnet")
library("IsingSampler")
library("parallel")
library("matrixStats")

source("helpers.R")
source("estimation.R")

# Set Simulation parameters -------------------------------------------
prob_interaction  <- .5    #c(.2, .4, .6)
N                 <- 250   #c(200, 1000, 5000)
p                 <- 5     #c(10, 25, 40)

iter              <- 1     #1000

precision_CI     <- .975  #c(.95, .997)
MCMC_draws        <- 1e5   #2.5e5

# Data Simulation -----------------------------------------------------
IM <- build_IM(p = p, prob = prob_interaction)  # build Ising model 

Data <- IsingSampler::IsingSampler(n = N, 
                                   graph = IM$Graph, 
                                   thresholds = IM$Thresholds)

# match true parameters to rbinnet model parameter output
true_parameters <- IM$Graph
diag(true_parameters) <- IM$Thresholds

# Estimation ----------------------------------------------------------
estimated_model <- estimate_IM(Data, 
                               MCMC_iter = MCMC_draws, 
                               precision = CI_precision)

# Performance ---------------------------------------------------------
corr <- cor(as.numeric(true_parameters), as.numeric(estimated_model$parameters))

# inclusion BF per edge
BF_inc <- estimated_model$inclusion_prob / (1 - estimated_model$inclusion_prob)

# sensitivity & specificity
SPE <- sum((BF_inc[lower.tri(BF_inc)] < 1) & (true_parameters[lower.tri(true_parameters)] == 0)) / sum(true_parameters[lower.tri(true_parameters)] == 0)
SEN <- sum((BF_inc[lower.tri(BF_inc)] > 1) & (true_parameters[lower.tri(true_parameters)] != 0)) / sum(true_parameters[lower.tri(true_parameters)] != 0)

# output measures
performance <- list(correlation = corr,
                    sensitivity = SEN,
                    specificity = SPE)

# Plots ---------------------------------------------------------------
plot(as.numeric(true_parameters), as.numeric(estimated_model$parameters))
