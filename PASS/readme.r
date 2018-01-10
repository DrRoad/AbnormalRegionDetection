#### This package contains 5 files. 
#### The PASS procedure is implemented in the function "pass.r".
#### Here is an example. 

rm(list = ls(all = TRUE));

source("pass.r");
source("hc.r");
source('peaks_new.R');

#### Input data
data = as.matrix(read.table('Data.txt'))

#### Set parameters in pass function. 
# Window size should be large enough to cover the longest segment
l=10
# Start point > 1 to stabilize the finite sample performance 
alpha0 = 5  
# The simulated threshold controls over-selection under null hypothesis
thres = 7

####  Call the pass function 
sel = pass(data, l, alpha0, thres); 

#### You will get a list of candidates for the recurrent variants in multiple samples
#### The results include the location of the segments and their corresponding test statistic 
#> sel
#       I   J     xstar
#[1,] 100 105 1332.9114
#[2,] 300 305 1084.4422
#[3,] 200 205  697.7702


#########################################

#### The data in Data.txt is generated as the follows. 
# 3 recurrent variants are generated at [100: 105], [200, 205], and [300:305]
# with carrier proportions 0.04, 0.06, and 0.08. 

N = 200
T = 500
mu = 2; prop= c(0.04, 0.06, 0.08)
q = length(prop);
s = rep(6, q);
loc = c(100, 200, 300);

data = matrix(0, N, T);
for (i in 1 : N){
    data[i,] = rnorm(T);
}
for (j in 1:q){
    for (i in 1:round(prop[j]*N)){
        data[i,loc[j]:(loc[j]+s[j]-1)] = data[i,loc[j]:(loc[j]+s[j]-1)]+mu;
    }
}

write.table(data, file = "Data.txt",  row.names = FALSE, col.names=FALSE)
