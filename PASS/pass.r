#### Identify recurrent segment variants from multiple samples
## A=data; L=l; a=alpha0; th=thres

pass = function(A, L, a, th) {
N = dim(A)[1]
T = dim(A)[2]

len = (T-L)*(L+1)+(L+1)*L/2;
x = rep(0, N);
HC = rep(0, len);

for (k in  1:len){
    t = ceiling(k/(L+1));
    j = min(k-(t-1)*L, T);
    for (i in 1:N) {
        x[i] = sum(A[i,t:j])/sqrt(j-t+1);
    }
    HC[k] = hc_2side(x, a)
} 
ind = which(HC>th);
can = HC[ind];
pe = peaks(ind, can, L);
select = pe[pe[,2]>pe[,1]+1,];
return(select); 
} 






