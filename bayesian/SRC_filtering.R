# calculates the filtering recursions and does the resampling
# we use a default value of alpha = 10^-4

SRC.filtering <- function( data , alpha = 1e-4 , affected_dim){ 

  # length and dimension of data
  n = dim(data)[1]
  N = dim(data)[2]
  # data summaries etc
  S <- rbind( rep(0,dim(data)[2]) , apply(data , 2 , cumsum) )
  S_2 <- cumsum( c( 0 , rowSums(data^2) ) )
  # p is used to calc the marginal like for abnormal
  # ratio of abnormal profiles
  p <- affected_dim/N

  # useful to define log of resampling probability
  log.alpha <- log(alpha)
  # lists for filtering, probs, location and types of segment
  weights <- vector("list",n)
  locations <- vector("list",n)
  type <- vector("list",n)
  # type - 0 for normal segment 1 for abnormal segment
  
  # initial N
  curr.locations <- c()
  curr.weights <- c()
  curr.type <- c()
  
  curr.locations[1] <- 0
  curr.weights[1] <- qN - 0.5 * ( S_2[2] - S_2[1] )
  curr.type[1] <- 0

  # initial A
  curr.locations[2] <- 0
  curr.weights[2] <- qA - 0.5 * ( S_2[2] - S_2[1] ) + P.A(1 , 1 , mu_seq)
  curr.type[2] <- 1

  c <- max(curr.weights)
  weights[[1]] <- curr.weights - ( c + log( sum( exp( curr.weights - c ) )  ) )
  locations[[1]] <- curr.locations
  type[[1]] <- curr.type
  
  for (t in 2:n){
  
    prev.type <- type[[t-1]]
    prev.locs <- locations[[t-1]]
    prev.weights <- weights[[t-1]]
    
    curr.type <- type[[t]]
    curr.locs <- locations[[t]]
    curr.weights <- weights[[t]]
   
    # label support points makes it easier to fnd positions in matrices
    sup.point.N <- which( prev.type == 0 )
    sup.point.A <- which( prev.type == 1 )
  
    # propagate normal particles
    # carry on in N state, nbinom.len
    # and P_N marginal likelihood
    
    i <- prev.locs[sup.point.N]
    nbinom.len <- pnbinom(t-i-1,k_N,p_N,log.p=T,lower.tail=F) - pnbinom(t-i-2,k_N,p_N,log.p=T,lower.tail=F) 
    curr.weights[sup.point.N] <- prev.weights[sup.point.N] + nbinom.len - 0.5 * ( S_2[ t+1 ] - S_2[ t ] )
    
    # propagate abnormal particles
    # carry on in A state, nbinom.len
    # Marginal likelihood PA calculated at i+1:t and i+1:t-1
    for (j in sup.point.A){
    
      i <- prev.locs[j]  
      nbinom.len <- pnbinom(t-i-1,k_A,p_A,log.p=T,lower.tail=F) - pnbinom(t-i-2,k_A,p_A,log.p=T,lower.tail=F) 
      curr.weights[j] <- prev.weights[j] + nbinom.len - 0.5 * ( S_2[ t+1 ] - S_2[ t ] ) + P.A( i+1 , t , mu_seq ) - P.A( i+1 , t-1 , mu_seq) 
    
    }
    
    ################################################
    # Calc other point chpt - C_t = t-1 , B_t = N/A
    ################################################
    
    ####
    # B_t = N
    ####
    # transition from abnormal to normal/ if there are no abnormal pts
    # give it a weight of zero i.e. log weight = - Inf
    i <- prev.locs[sup.point.A] 

    if( length(i)==0 ){
      C_N <- - Inf
    }
    else{ 
      # weight to propogate
      W <- prev.weights[sup.point.A] + dnbinom(t-i-1,k_N,p_N,log=T) - pnbinom(t-i-2,k_N,p_N,log.p=T,lower.tail=F) + log(pi.N)
      
      # transition from A - N
      # (only trans possible to get to N) 
      # prob of pi.N
      cmax <- max( W )
      C_N <- - 0.5*( S_2[ t+1 ] - S_2[ t  ] ) + cmax + log( sum( exp(W - cmax) ) ) 
    }
  
    ####
    # B_t = A
    ####
    # transition to abnormal, can be A-A or N-A
    i.N <- prev.locs[sup.point.N]
    i.A <- prev.locs[sup.point.A]
    
    if ( length(i.N) == 0 ){
      
      # no normal particles just use A particles
      W <- prev.weights[sup.point.A] + dnbinom( t - i.A -1  ,k_A,p_A,log=T) - pnbinom(t - i.A - 2,k_A,p_A,log.p=T,lower.tail=F) 
      cmax <- max( W )
      C_A <- P.A(t,t,mu_seq) - 0.5*( S_2[ (t+1) ] - S_2[ t  ] ) + log(pi.A) + cmax + log( sum( exp(W - cmax) ) ) 
    
    }
    
    else if ( length(i.A) == 0 ) {
      
      # no abnormal particles just use N particles
      W <- prev.weights[sup.point.N] + dnbinom( t - i.N -1 ,k_N,p_N,log=T) - pnbinom(t - i.N - 2,k_N,p_N,log.p=T,lower.tail=F) 
      cmax <- max( W )
      C_A <- P.A(t,t,mu_seq) - 0.5*( S_2[ (t+1) ] - S_2[ t  ] ) + cmax + log( sum( exp(W - cmax) ) ) 
      
    } 
    
    else {
      
      # each have some N and A particles do seperately
      W <- prev.weights[sup.point.A] + dnbinom( t - i.A-1  ,k_A,p_A,log=T) - pnbinom(t - i.A - 2 ,k_A,p_A,log.p=T,lower.tail=F) 
      cmax <- max( W )
      ABS.PART <- log(pi.A) + cmax + log( sum( exp(W - cmax) ) ) 
      W <- prev.weights[sup.point.N] + dnbinom( t - i.N -1  ,k_N,p_N,log=T) - pnbinom(t - i.N - 2 ,k_N,p_N,log.p=T,lower.tail=F) 
      cmax <- max( W )
      NORM.PART <- cmax + log( sum( exp(W - cmax) ) ) 
      # put these together 
      part <- c( NORM.PART , ABS.PART )
      cmax <- max(part)
      C_A <- P.A(t,t,mu_seq) - 0.5*( S_2[ (t+1) ] - S_2[ t  ] ) + cmax + log( sum( exp( part - cmax ) ) )
      
    }
 
    ###############################################
     ## resample if neccessary only when t > 20 ##
    ###############################################
  
    if ( t > 20 ){
    
      temp.weights <- c( curr.weights , C_N , C_A )
      temp.locs <- c( prev.locs , (t-1) , (t-1) )
      temp.type <- c( prev.type ,  0 , 1 )
    
      ## log normalized weights ##
      c <- max( temp.weights )
      lognorm.weights <-  temp.weights - ( c + log( sum(exp(temp.weights-c)) ) ) 
    
      ##############################
      # stratified resampling part #
      ##############################
      # take all the log weights that are < alpha and resample them
      to.resample <- lognorm.weights[ lognorm.weights < log.alpha ]
      # resampling function in resample.r
      lognorm.weights[ lognorm.weights < log.alpha ] <-  resample( to.resample , alpha )
    
      # remove weights for which = -Inf
      lognorm.weights[lognorm.weights == -Inf] <- NA
      updated.locs <- temp.locs[ !is.na(lognorm.weights) ]
      updated.type <- temp.type[ !is.na(lognorm.weights) ]
    
      # renormalise weights - is this necessary?
      temp.weights <- lognorm.weights[!is.na(lognorm.weights)]
      c <- max( temp.weights )
      renorm.weights <-  temp.weights - ( c + log( sum(exp(temp.weights-c)) ) ) 
    
      # put back in vectors (ordered)
      weights[[t]] <- renorm.weights 
      locations[[t]] <- updated.locs 
      type[[t]] <- updated.type 
    
    }
  
    #########################################
         ## No resampling t <= 20 ##
    #########################################
  
    else{
      
      temp.weights <- c( curr.weights , C_N , C_A )
      
      # find log normalised weights #
      c <- max( temp.weights )
      lognorm.weights <-  temp.weights - ( c + log( sum(exp(temp.weights-c)) ) )
      
      weights[[t]] <- lognorm.weights 
      locations[[t]] <- c( prev.locs , (t-1) , (t-1) )
      type[[t]] <- c( prev.type , 0 , 1  )
    
    }
  
  }
  
  newList <- list("weights" = weights , "locations" = locations , "type" = type )
  return(newList)

}

