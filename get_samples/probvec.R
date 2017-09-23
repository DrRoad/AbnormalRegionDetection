get.probvec <- function( filtering , no.draws=1000 ){

  logprobs <- filtering$weights
  cpt.locs <- filtering$locations
  types <- filtering$type

  # vec is a vector of counts of seg locations #
  vec <- numeric(n)
  prob.states <- matrix(nrow=no.draws,ncol=n)
  
  for (i in 1:no.draws){
  
    f <- draw.from.post(logprobs,cpt.locs,types)
    segs <- f$draws
    states <- f$states
  
    # fill vec (hist of chpts)
    vec[f$draws] <- vec[f$draws] + 1
  
    # probs of being N or A
    prob.states[i,] <- get.states( segs , states )
    
  }

  return( apply(prob.states,2,sum)/no.draws )

}