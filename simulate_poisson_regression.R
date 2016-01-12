ndot=500
cps=c(.2*ndot,.3*ndot,.7*ndot)
dcps=diff(c(0,cps,ndot))
means=c(1,3,5,2)
x<-NULL
for(i in 1:length(dcps)) x<-c(x,rpois(dcps[i],means[i]))
write(x,'pr.txt',ncolumn=1)#length(x))