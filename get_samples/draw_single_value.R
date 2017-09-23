# draw one sample of chpts and states from the posterior

draw.from.post <- function( logprobs , cpt.locs , types ){ 

  n <- length(logprobs)
  l.probs <- logprobs[[n]]
  seg.locs <- cpt.locs[[n]]
  seg.types <- types[[n]]

  ind <- 1:length(l.probs)
  a <- sample( ind , size=1 , replace = TRUE, exp( l.probs ) )

  # state amd location
  curr.state <- seg.types[a]
  t <- seg.locs[a]

  STATES <- curr.state
  DRAW <- t
  
  while (t > 1){

    if(curr.state==0){
    
      l.probs <- logprobs[[t]][ types[[t]] == 1 ] 
      i <- cpt.locs[[t]][ types[[t]] ==1 ]
      back.dens <- log(pi.N) + dnbinom(t-i-1,k_A,p_A,log=TRUE) - pnbinom(t-i-2,k_A,p_A, lower.tail=FALSE , log.p=TRUE)
      l.probs <- l.probs + back.dens  
    
      ind <- 1:length(l.probs)
      a <- sample( ind , size=1 , replace = TRUE, exp( l.probs ) )
    
      # state amd location
      curr.state <- 1
      t <- i[a]
      STATES <- c(STATES,curr.state)
      DRAW <- c(DRAW,t)
    }
    else{
    
      # currently in abnormal state poss transitions to normal or abnormal
      i.N <- cpt.locs[[t]][ types[[t]] == 0 ]
      i.A <- cpt.locs[[t]][ types[[t]] == 1 ]
    
      l.probsA <- logprobs[[t]][  types[[t]] == 1 ] 
      back.densA <- log(pi.A) + dnbinom(t-i.A-1,k_A,p_A,log=TRUE) - pnbinom(t-i.A-2,k_A,p_A, lower.tail=FALSE , log.p=TRUE)
      l.probsA <- l.probsA + back.densA
      # normal points
      l.probsN <- logprobs[[t]][  types[[t]] == 0 ] 
      back.densN <- dnbinom( t-i.N-1 , k_N , p_N , log=TRUE ) - pnbinom( t-i.N-2 , k_N , p_N , lower.tail=FALSE , log.p=TRUE )
      l.probsN <- l.probsN + back.densN
      # now normalise
      i.all <- c( i.N , i.A )
      t.all <- c( rep(0,length(i.N)) , rep(1,length(i.A)) )
      l.probs <- c( l.probsN , l.probsA )
      c <- max(l.probs)
      l.probs.norm <- l.probs - ( c + log( sum( exp( l.probs - c ) ) ) ) 
    
      ind <- 1:length(l.probs.norm)
      a <- sample( ind , size=1 , replace = TRUE, exp( l.probs.norm ) )
    
      # state amd location
      curr.state <- t.all[a]
      t <- i.all[a]
      STATES <- c(STATES,curr.state)
      DRAW <- c(DRAW,t)
      }
  
   }
  
  DRAW[DRAW==0] <- 1
  newlist <- list("draws"=c(n,DRAW),"states"=STATES)
  return(newlist)

}
  