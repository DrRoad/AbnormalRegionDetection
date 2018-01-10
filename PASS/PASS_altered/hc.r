# Calcualte the HC statistic based on a vector of standardized measures and a starting point 
# Assume normal mixure model

####### 2-sided HC
hc_2side = function(y, a) {
N = length(y)
p_val = 2*pnorm(abs(y), lower.tail=F) # 2-sided test
sortp = sort(p_val);
ind = c(a: (N/2)) 
tmp = (ind/N - sortp[ind])/sqrt(sortp[ind]*(1-sortp[ind]))
peak= sqrt(N)*max(tmp)
return(peak) 
} 

###### 1-sidedm (right) HC
hc_right = function(y, a) {
N = length(y)
p_val = pnorm(y, lower.tail=F) # right-sided test
sortp = sort(p_val);
ind = c(a: (N/2)) 
tmp = (ind/N - sortp[ind])/sqrt(sortp[ind]*(1-sortp[ind]))
peak= sqrt(N)*max(tmp)
return(peak) 
}

###### 1-sidedm (left) HC
# y = b[,x]; a=start
hc_left = function(y, a) {
N = length(y)
p_val = pnorm(y, lower.tail=T) # right-sided test
sortp = sort(p_val);
ind = c(a: (N/2)) 
tmp = (ind/N - sortp[ind])/sqrt(sortp[ind]*(1-sortp[ind]))
peak= sqrt(N)*max(tmp)
return(peak) 
}
