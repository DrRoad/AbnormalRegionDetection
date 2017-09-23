# apply loss function to get a segmentation returns vector with 0's and 1's
# 0 means in a normal segment , 1 is abnormal
# apply loss with parameter gamma

loss <- function( gamma, probvec ){

  p.gamma <- 1/(1+gamma)
  probvec <- prob.vec
  probvec[probvec < p.gamma ] <- 0
  probvec[probvec>0] <- 1
  return(probvec)
  
}
  
