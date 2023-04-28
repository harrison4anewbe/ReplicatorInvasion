source("invasion_graph_main_functions.R")
source('extrafunctions.R')
library("deSolve")

Replicator.IS=function(P,tolerance=1e-14){
  k=dim(P)[1] # number of species
  C=list() # list to hold all the feasible (i.e. non-negative entries) equilibria of the model
  for (i in 1:k){ # add first k equilibrea where the are of the form (1,0,0,..0), (0,1,0,...,0) etc
    temp = rep(0,k)
    temp[i] = 1
    C = append(C, list(temp))
  }
  no.C=k # counter to keep track of the number of feasible equilibria
  # find all the feasible equilibria
  for(i in 2:k){ #loop for the number of species supported by the equilibrium. Note: for replicator this will start with 2 instead of 1
    temp=combn(1:k,i) # a matrix of all the possible configurations of communities with i species
    k2=dim(temp)[2] # the number of configurations
    for(j in 1:k2){ # loop through the configurations
      I=temp[,j] # pull out the j-th configuration
      if (abs(det(P[I,I]))<tolerance){ # check if determinant of submatrix is 0
        tempor = which(abs(eigen(P[I,I])$values)<tolerance) # identify which eigenvalues are 0
        for(x in tempor){
          vec = Re(eigen(P[I,I])$vectors[,x])/sum(Re(eigen(P[I,I])$vectors[,x])) # normalize eigenvector corresponding to the eigenvalue of 0
          boolvec <- sign(Re(vec)) == rep(1,dim(P[I,I])[1]) # vector of TRUEs and FALSEs telling if all elements of eigenvector are same sign
          if (identical(boolvec, rep(TRUE, dim(P[I,I])[1]))==TRUE){ # check if all elements of eigenvector are same sign
            xstar = vec # let xstar be the normalized eigenvector
            no.C=no.C+1 # update counter
            xstar2=rep(0,k) # grab the zero vector to fill out the feasible equilibrium for all k species
            xstar2[I]=xstar # fill in the non-zero entries
            C[[no.C]]=xstar2 # set as the next element in the list
          }
        }
      } else{
        xtemp=solve(a=P[I,I],b=rep(1,length(I))) # solve for the equilibrium restricted to the species in I
        xstar = xtemp/sum(xtemp)
        if(min(c=xstar)>0){ # make sure all entries are positive i.e. a new feasible equilibrium
          no.C=no.C+1 # update counter
          xstar2=rep(0,k) # grab the zero vector to fill out the feasible equilibrium for all k species
          xstar2[I]=xstar # fill in the non-zero entries
          C[[no.C]]=xstar2 # set as the next element in the list
          stop
        }
      }
    }
  }
  # create the invasion scheme matrix
  IS=matrix(NA,length(C),k) # matrix to hold all the per-capita growth rates
  for(i in 1:length(C)){
    IS[i,]=P%*%C[[i]]-(t(C[[i]])%*%P%*%C[[i]])[1,1] # i-th row corresponds to all the per-capita growth rates at the i-th feasible equilibrium
  }
  IS[which(abs(IS)<tolerance)]=0 # set to zero based on the tolerance
  return(IS)
}
