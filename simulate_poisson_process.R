set.seed(0)
end=1e4
ncp=500
cps=c(0,sort(runif(ncp))*end)
dcps=diff(c(cps,end))
lambdas=rgamma(ncp+1,1,0.01)
ndata=rpois(ncp+1,lambdas*dcps)
x=NULL
for(i in 1:(ncp+1)) x=c(x,cps[i]+dcps[i]*sort(runif(ndata[i])))
write(x,'pp.txt',ncolumn=1)