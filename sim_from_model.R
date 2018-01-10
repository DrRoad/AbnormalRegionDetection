# simulate from the model -- data like from the first simulation study 
# in the paper, dimension 200 and length approx 1000

# n.max is a soft maximum on length of data set 
n.max <- 1000

# fixed no of dimensions and no of affected dimensions i.e. 4% 
N = 200
affected_dim <- 8

# start in a normal state encoded as 0 
cpts <- rnbinom( 1 , k_N , p_N )
state <- rep( "N", cpts )

data <- matrix( rnorm(cpts*N) , nrow = cpts ,ncol = N)
curr.state <- 0 
n <- cpts

while (n < n.max){

  # if current state is normal - draw abnormal length
  if (curr.state == 0){
    
    # in normal transition to abnormal length ~ NBinom( k_A , p_A )
    seg.length <- rnbinom( 1 , k_A , p_A )
    cpts <- c( cpts , tail(cpts,1)+seg.length )
    state <- c( state , rep( "A" , seg.length ) )
    
    # grow data seg.length no of rows
    data <- rbind( data, matrix(rnorm(seg.length*N),ncol=N) )
    
    # draw value of mu for dims with affected means
    mu <- runif(1 , 0.3 , 0.7) 
    data[ (tail(cpts,2)[1]+1):tail(cpts,1) ,  1:affected_dim  ] <- data[ (tail(cpts,2)[1]+1):tail(cpts,1) , 1:affected_dim ] + mu 
    
    # current state abnormal update it
    curr.state <- 1
    
  }
  
  else{
    
    # current state is abnormal curr.state == 1
    # with prob pi.N go to Normal
    # with prob pi.A go to Abnormal with diff mu
    # generate a u bernoulli with probs pi.A,pi.N
    u <- rbinom(1,1,pi.A) 
    
    # u==0 transition to normal
    
    if (u==0){
      
      seg.length <- rnbinom( 1 , k_N , p_N )
      cpts <- c( cpts , tail(cpts,1)+seg.length )
      state <- c( state , rep( "N" , seg.length ) )
      
      # grow data seg.length no of rows
      data <- rbind( data, matrix(rnorm(seg.length*N),ncol=N) )
      
      curr.state <- 0
      
    }
    else{
      
      # u==1 thus transition to abnormal
      seg.length <- rnbinom( 1 , k_A , p_A )
      cpts <- c( cpts , tail(cpts,1)+seg.length )
      state <- c( state , rep( "A" , seg.length ) )
      
      # grow data seg.length no of rows
      data <- rbind( data, matrix(rnorm(seg.length*N),ncol=N) )
      
      # draw mu from prior
      mu <- runif(1 , 0.3 , 0.7) 
      
      # random dimensions to affect ran.affec
      # ran.affec is first one
      ran.affec <- sample(affected_dim:(N-affected_dim),1)
      data[ (tail(cpts,2)[1]+1):tail(cpts,1) ,  ran.affec:(ran.affec+affected_dim)  ] <- data[ (tail(cpts,2)[1]+1):tail(cpts,1) , ran.affec:(ran.affec+affected_dim) ] + mu 
      
      curr.state <- 1
      
    }
 
    
  }  
  
  n <- tail(cpts,1)
  
}


## stat distribution qN , qA
EN <- ( k_N * (1-p_N) )/p_N
EA <- ( k_A * (1-p_A) )/p_A

ldenom <- log( pi.N*EN + EA )
qA <- log(pi.N) + log(EN) - ldenom 
qN <- log(EA) - ldenom

# prior for mu 
# gridded to calculate 
# marginal like for abnormal P_A in PAC(.,)
mu_seq <- seq(0.3 , 0.7 , by=0.05 )

S <- rbind( rep(0,dim(data)[2]) , apply(data , 2 , cumsum) )
S_2 <- cumsum( c( 0 , rowSums(data^2) ) )
# p is used to calc the marginal like for abnormal
# ratio of abnormal profiles
p <- affected_dim/N
