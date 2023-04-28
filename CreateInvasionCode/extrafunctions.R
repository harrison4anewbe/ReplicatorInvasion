# functions to (1) find cycles and (2) to understand what happens at the end of community assembly

######################
# THE FindCycles FUNCTION
######################
# Command for finding cycles in a cyclic graph
# Source: https://stackoverflow.com/a/55094319
# Input: a graph (igraph format)
# Output: some (but not necessarily all) of the directed cycles in the graph (see https://stackoverflow.com/a/55094319 for more details)
FindCycles = function(g) {
  Cycles = NULL
  for(v1 in V(g)) {
    if(degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
      TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  Cycles
}

########
# end states - assuming acyclic
# input is the output of IG.function
########

end.state=function(out){
  k=dim(out$IS)[2] # number of species
  # go through all minus i communities and save if an end state
  end.state=c()
  for(i in 1:k){
    for(j in out$minus.i[[i]]){
      if(out$IS[j,i]<0)end.state=c(end.state,j)
    }
  }
  if(length(end.state)==0){end.state=dim(out$IS)[1]}
  return(end.state)
}

# GLV ODE solver and plotter

GLV_plotter=function(A,b,initialx,times=1:100)
{
  lv=function(t,x,parms){
    with(parms,{
      dx<-x*(b+A%*%x)
      list(dx)
    })
  }
  k=length(b)
  parms=list(b=b,A=A)
  out=ode(y=initialx,times=times,func=lv,parms=parms)
  par(mar=c(4.5,4.5,0,1))
  matplot(out[,1],out[,-1],type="l",lty=1,lwd=4,bty="n",xlab="time t",ylab=expression(paste("densities ",x[i])))
}



# LV SDE simulator
GLV_SDE_plotter=function(A,b,initialx,Tend=100,h=0.05,sigma=0.25,reps=1){
  n=length(b)
  steps=floor(Tend/h)
  x=array(NA,dim = c(steps+1,n,reps))
  for(j in 1:reps){
  x[1,,j]=initialx
  for(i in 1:steps){
    U=(A%*%x[i,,j]+b-sigma^2/2)*h+sigma*sqrt(h)*rnorm(n)
    x[i+1,,j]=x[i,,j]*exp(U)
  }
  }
  par(mar=c(4.5,4.5,0,1))
  times=seq(0,Tend,by=h)
  matplot(times,x[,,1],type="l",lty=1,lwd=3,bty="n",xlab="time t",ylab=expression(paste("densities ",X[i])))
}
  


