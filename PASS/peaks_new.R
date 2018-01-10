# Locate signal segments.

peaks = function(index, stat, bin){

left = ceiling(index/(bin+1));
right = index-(left-1)*bin;

list = cbind(left, right, stat)

I = rep(0,1);
J = rep(0,1);
xstar = rep(0,1);

t=1;

while (length(list)> 3) {

k.max = which(abs(list[,3])== max(abs(list[,3])));
I[t] = list[k.max,1]
J[t] = list[k.max,2]
xstar[t] = list[k.max,3];

II = which(list[,1]<=J[t] & list[,2]>=I[t]);

list = list[-II,];

t=t+1;
}

if (length(list)==3){
I[t] = list[1]; 
J[t] = list[2]
xstar[t] = list[3];
}

location = cbind(I, J, xstar);
return(location);
}
