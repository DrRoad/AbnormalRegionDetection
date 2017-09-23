##########################
# Example use of BARD
###########################
rm(list=ls())
# to make a little faster
#require(compiler)
#enableJIT(3)

###########################
# filtering functions
###########################
# function to compute abnormal marginal likelihood
source("bayesian/PA.R")
# resampling function for SRC
source("bayesian/resample.R")
# full filtering and SRC function
source("bayesian/SRC_filtering.R")
# cutoff for SRC filtering 
alpha <- 1e-4

###########################
# Sampling from the posterior and loss 
###########################
source("get_samples/draw_single_value.R")
source("get_samples/get_states.R")
source("get_samples/probvec.R")
source("get_samples/loss.R")

############################
# simulate from model
############################
source("params.R")
source("sim_from_model.R")

# calculate posterior recursions with resampling
filtering <- SRC.filtering(data,alpha,affected_dim) 

# draw from the posterior no.draws times
no.draws = 1000
prob.vec <- get.probvec(filtering,no.draws)
# returns a vector of length data with prob.vec[t] as the probability of time point t
# being abnormal

# this can be plotted etc
plot(prob.vec)
# add lines for real chpts
abline(v=cpts)

# to get an explicit segmentation - 0's for Normal , 1's for Abnormal
gamma = 1/3
segmentation <- loss(gamma,prob.vec)

