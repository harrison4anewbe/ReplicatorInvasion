# Code to plot zero sum replicator equation
# A is a k x k payoff matrix where k is the number of species

library(deSolve)

RE_plotter=function(A,initialx,times=1:100)
{
  re=function(t,x,parms){
    with(parms,{
      dx<-x*(A%*%x)
      list(dx)
    })
  }
  k=length(A[1,])
  parms=list(A=A)
  out=ode(y=initialx,times=times,func=re,parms=parms)
  par(mar=c(4.5,4.5,0,1))
  matplot(out[,1],out[,-1],type="l",lty=1,lwd=4,bty="n",xlab="time t",ylab=expression(paste("densities ",x[i])), col = cols)
  #legend(1, 12, spp, col = cols, pch =15)
}