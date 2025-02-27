for(i in 1:42){cat(paste("set.seed(1234)

library(igraph)
library(tidyverse)
library(mgcViz)
library(ggplot2)
library(parallel)
library(gridExtra)
library(scales)
library(coin)
library(expandFunctions)
library(grid)



get.B0hat <- function(A, K, block_memberships) {
  B0hat = matrix(0, K, K)
  n.lk = matrix(0, K, K)
  for (l in 1:K) {
    for (k in 1:l) {
      temp = A[block_memberships == l, block_memberships == k]
      edges = sum(temp)
      size = length(c(temp))
      if (l == k) {
        n.lk[l, k] = sum(upper.tri(temp, diag = F))
        B0hat[l, k] = edges / (2 * n.lk[l, k])
      } else{
        n.lk[l, k] = n.lk[k, l] = size
        B0hat[l, k] = B0hat[k, l] = edges / size
      }
    }
  }
  return(list(B0hat, n.lk))
}


get.B1hat <- function(A, K, block_memberships, stencil, n.lk) {
  Nlk.diag <- rep(0, K)
  B1.diag <- rep(0, K)
  
  for (l in unique(diag(stencil))) {
    blocks <- which(stencil == l, arr.ind = T)[, 1]
    edges = 0
    size = 0
    for (b in blocks) {
      temp = A[block_memberships == b, block_memberships == b]
      edges = edges + sum(temp) / 2
      size = size + n.lk[b, b]
    }
    Nlk.diag[blocks] = size
    B1.diag[blocks] = edges / size
  }
  N.lk = matrix(0, K, K)
  B1 = matrix(0, K, K)
  for (l in unique(stencil[lower.tri(stencil)])) {
    blocks.1 = which(stencil == l, arr.ind = T)[, 1]
    blocks.2 = which(stencil == l, arr.ind = T)[, 2]
    edges = 0
    size = 0
    for (b in 1:length(blocks.1)) {
      temp = A[block_memberships == blocks.1[b], block_memberships == blocks.2[b]]
      edges = edges + sum(temp)
      size = size + length(c(temp))
      B1[blocks.1[b], blocks.2[b]] = B1[blocks.2[b], blocks.1[b]] = NA
    }
    N.lk[which(is.na(B1))] = size
    B1[which(is.na(B1))] = edges / size
  }
  diag(N.lk) = Nlk.diag
  diag(B1) = B1.diag
  return(list(B1, N.lk))
}


LLRind <- function(B0hat, B1hat, n.lk, stencil) {
  K=dim(n.lk)[1]
  temp=n.lk * B0hat * log(B0hat / (1 - B0hat)) +
    n.lk * log(1 - B0hat) -
    n.lk * B0hat * log(B1hat / (1 - B1hat)) -
    n.lk * log(1 - B1hat)
  temp[which(is.nan(temp), arr.ind = T)]=0
  
  return(temp)
}


gamma.B=function(stencil){
  K=dim(stencil)[1]
  gammaB=matrix(0,K,K)
  for( l in 1:K){
    for(k in 1:l){
      if(l==k){
        gammaB[l,k]=sum(stencil==stencil[l,k])-1
      }else{
        gammaB[l,k]=sum(stencil==stencil[l,k])/2-1
      }
      
    }
  }
  gammaB[upper.tri(gammaB)]=t(gammaB)[upper.tri(gammaB)]
  return(gammaB)
  
}

get.stencil.general=function(Tree.full,nD=0,nZ=0){
  stencils.inner=list()
  skip=0
  param.offset=0
  for(i in Tree.full$clct){
    Tree=list(motif.sizes=Tree.full$motif.sizes,
              motifs=Tree.full$motifs[(skip+1):(skip+i)])
    skip=skip+i
    stencils.inner=append(stencils.inner,list(get.stencil.inner(Tree,param.offset)))
  }
  
  total.param=max(unlist(stencils.inner))
  
  stencil=c()
  
  for(i in 1:length(Tree.full$clct)){
    M=c()
    for(j in 1:length(Tree.full$clct)){
      if(i>j){
        dmnts=c(dim(stencils.inner[[j]])[1],dim(stencils.inner[[i]])[2])
        M=rbind(matrix(0,dmnts[1],dmnts[2]),M)
      }else if(i==j){
        M=rbind(M,stencils.inner[[i]])
      }else if(i<j){
        total.param=total.param+1
        dmnts=c(dim(stencils.inner[[j]])[1],dim(stencils.inner[[i]])[2])
        M=rbind(M,matrix(total.param,dmnts[1],dmnts[2]))
      }
    }
    stencil=cbind(stencil,M)
  }
  
  stencil[upper.tri(stencil)]= t(stencil)[upper.tri(stencil)]
  total.param=max(stencil)
  
  if(nD==0) return(stencil)
  pD=sample(1:total.param,nD)
  for(i in pD){
    total.param=total.param+1
    reps=floor(length(which(stencil[lower.tri(stencil,diag = T)]==i))/10)
    sper=max(1,reps)
    ind=sample(which(stencil[lower.tri(stencil,diag = T)]==i),sper)
    stencil[lower.tri(stencil, diag = T)][ind]=total.param
  }
  stencil[upper.tri(stencil)]= t(stencil)[upper.tri(stencil)]
  
  return(list(stencil,pD))
}

get.stencil.inner=function(Tree,param.offset=0){
  motif.stencil=list()
  
  total.param=0
  
  for(comm.size in Tree$motif.sizes){
    M=matrix(0,comm.size,comm.size)
    M[lower.tri(M,diag = T)]=1:(choose(comm.size, 2) + comm.size)+total.param
    total.param=max(M)
    motif.stencil=append(motif.stencil,list(M=M))
  }
  init.param=total.param
  total.param=total.param+param.offset
  stencil=c()
  for(i in 1:length(Tree$motifs)){
    M=c()
    for(j in 1:length(Tree$motifs)){
      if(i>j){
        dmnts=Tree$motif.sizes[Tree$motifs[c(i,j)]]
        M=rbind(matrix(0,dmnts[2],dmnts[1]),M)
      }else if(i==j){
        M=rbind(M,motif.stencil[[Tree$motifs[i]]])
      }else if(i<j){
        total.param=total.param+1
        dmnts=Tree$motif.sizes[Tree$motifs[c(i,j)]]
        M=rbind(M,matrix(total.param,dmnts[2],dmnts[1]))
      }
    }
    stencil=cbind(stencil,M)
  }
  
  stencil[upper.tri(stencil)]= t(stencil)[upper.tri(stencil)]
  
  param.offset<<-total.param-init.param
  
  return(stencil)
}

get.B=function(stencil,alphas){
  K=dim(stencil)[1]
  nparam=max(stencil)
  B=matrix(0,K,K)
  x=sample_dirichlet(nparam,alphas)
  for(i in 1:nparam){
    B[which(stencil==i,arr.ind = T)]=sum(x[,i]*x[,i])
    
  }
  return(B)
}

make.graph = function(n, B, block.sizes) {
  graph_sbm = sample_sbm(n,
                         B,
                         block.sizes = block.sizes,
                         directed = FALSE,
                         loops = FALSE)
  A = matrix(as_adjacency_matrix(graph_sbm), nrow = n)
  return(A)
}

get.t = function(vec) {
  temp = c()
  for (i in 1:length(vec)) {
    temp = c(temp, rep(i, vec[i]))
  }
  return(temp)
}


LLR <- function(B0hat, B1hat, n.lk) {
  L = dim(B0hat)[1]
  stat = 0
  counter = 1
  for (l in 1:L) {
    for (k in 1:l) {
      B0 = B0hat[l, k]
      B1 = B1hat[l, k]
      n = n.lk[l, k]
      if (B0 == 0 || B1 == 0) {
        
      } else{
        stat = stat + n * B0 * log(B0 / (1 - B0)) +
          n * log(1 - B0) -
          n * B0 * log(B1 / (1 - B1)) -
          n * log(1 - B1)
        counter = +1
        
      }
      
    }
  }
  return(stat)
}

t.error=function(rejected, pD, stencil){
  t1=0
  t2=0
  for(i in 1:max(stencil)){
    if(i %in% pD){
      t2=t2+rejected[which(stencil==i)][1]
    }else{
      t1=t1+rejected[which(stencil==i)][1]
    }
    
  }
  t2=1-t2/length(pD)
  t1=t1/(max(stencil)-length(pD))
  return(c(t1,t2))
}

BIC <- function(B0, B1, stencil, n.lk, N.lk) {
  penalty = 0
  for (i in 1:max(stencil)) {
    if (i != 0) {
      penalty = penalty + (sum(stencil == i) - 1)
    }
  }
  return(penalty)
}

anova.test.lapply=function(x){
  return(summary(aov(P~site,x))[[1]][['Pr(>F)']][1])
}

BH.rejection=function(P,stencil){
  P.ordered=P[order(P)]
  
  param=1:length(P)
  param.ordered=param[order(P)]
  regline=(1:length(P))*rep(0.05/length(P))
  k=suppressWarnings(max(which(P.ordered<regline)))
  if(k!=-Inf){
    ind=which(stencil %in% param.ordered[1:k])
  }else{
    ind=c()
  }
  
  rejected=matrix(0,dim(stencil)[1],dim(stencil)[1])
  rejected[ind]=1
  rejected[upper.tri(stencil)]=t(rejected)[upper.tri(stencil)]
  
  X=data.frame(P0=P.ordered,n=1:length(P),r=regline)
  
  return(list('rejected'=rejected, 'X'=X))
}

friedman.test.lapply.lr.mult= function(x){
  return(friedman.test(x$mat)$p.value)
}



m=",i,"

sims=readRDS(file='models.RData')

MODEL=sims[[m]]$model
NDS=sims[[m]]$nD
Tree.full=MODEL

FILESAVENAME=paste('model-',MODEL$name,'nDs-',NDS,'.RData', sep = '')

org.stencil=sims[[m]]$org.stencil
pD=sims[[m]]$pD

sample.outer=200
sample.inner=20
n.star=200




results.list.outer=list()


for( i in 1:sample.outer){
  results.list.inner=list()
  stencil=get.stencil.general(MODEL,nD=0)
  
  K=sum(MODEL$motif.sizes[MODEL$motifs])
  df = choose(K, 2) + K - max(stencil)
  stokes = qchisq(0.95, df)
  
  B=get.B(org.stencil,rep(1/dim(org.stencil)[1],dim(org.stencil)[1]))
  
  
  
  for(j in 1:sample.inner){
    corruption=matrix(rnorm(K^2,0,0.1),K,K)
    corruption=corruption*t(corruption)
    fix.diag=rbinom(length(diag(stencil)),1,prob=0.5)
    fix.diag[which(fix.diag==0)]=-1
    diag(corruption)=fix.diag*diag(corruption)
    
    block.sizes=rep(n.star,dim(B)[1])
    
    A=make.graph(sum(block.sizes),pmin(B+B*corruption,0.99),block.sizes)
    
    block.membership=get.t(block.sizes)
    
    B0hat.nlk=get.B0hat(A,length(block.sizes),block.membership)
    
    B0hat=B0hat.nlk[[1]]
    n.lk=B0hat.nlk[[2]]
    
    
    
    B1hat.Nlk=get.B1hat(A,length(block.sizes),block.membership,stencil,n.lk)
    
    B1hat=B1hat.Nlk[[1]]
    N.lk=B1hat.Nlk[[2]]
    
    gammaB=gamma.B(stencil)
    stat=2*LLRind(B0hat,B1hat,n.lk,stencil)
    P=pchisq(stat,gammaB,lower.tail = FALSE)
    
    P0=c()
    
    for(s in 1:max(stencil)){
      P0=c(P0,P[which(stencil==s)][1])
    }
    P0.ordered.chisq=P0[order(P0)]
    param=1:max(stencil)
    param.ordered=param[order(P0)]
    regline=(1:length(P0))*rep(0.05/length(P0))
    k=suppressWarnings(max(which(P0.ordered.chisq<regline)))
    if(k!=-Inf){
      ind=which(stencil %in% param.ordered[1:k])
    }else{
      ind=c()
    }
    
    G=1:dim(B)[1]
    
    rejectedBH.chisq=matrix(0,length(G),length(G))
    rejectedBH.chisq[ind]=1
    rejectedBH.chisq[upper.tri(P)]=t(rejectedBH.chisq)[upper.tri(P)]
    
    
    bonferroni=qchisq(0.05/(max(stencil)),gammaB)
    rejectedBF=matrix(0,length(G),length(G))
    rejectedBF[which(P<bonferroni)]=1
    rejectedBF[upper.tri(P)]=t(rejectedBF)[upper.tri(P)]
    
    llr = LLR(B0hat, B1hat, n.lk)
    
    t1.t2=t.error(rejectedBH.chisq,pD,stencil)
    
    t1.rejectedBH.chisq=t1.t2[1]
    t2.rejectedBH.chisq=t1.t2[2]
    
    t1.t2=t.error(rejectedBF,pD,stencil)
    
    t1.rejectedBF=t1.t2[1]
    t2.rejectedBF=t1.t2[2]
    
    n.s=sum(MODEL$motif.sizes[MODEL$motifs])*n.star
    llr.bic=llr-BIC(B0hat, B1hat, stencil, n.lk, N.lk)*log(choose(n.s,2))
    results.list.inner=append(results.list.inner,
                        list(list(Tree=MODEL$name,
                                  nD=NDS,
                                  rejectedBH.chisq=rejectedBH.chisq,
                                  rejectedBF.chisq=rejectedBF,
                                  llr=llr,
                                  llr.test=2*llr<stokes,
                                  llr.bic=llr.bic,
                                  B0=B0hat,
                                  B1=B1hat,
                                  P0=P0.ordered.chisq,
                                  n.lk=n.lk,
                                  N.lk=N.lk,
                                  t1.rejectedBH.chisq=t1.rejectedBH.chisq,
                                  t2.rejectedBH.chisq=t2.rejectedBH.chisq,
                                  t1.rejectedBF=t1.rejectedBF,
                                  t2.rejectedBF=t2.rejectedBF)))
    
  }
  
  stencil=get.stencil.general(Tree.full,nD=0)
  
  rejectedBH.avg=matrix(0,K,K)
  rejectedBF.avg=matrix(0,K,K)
  B0=matrix(0,K,K)
  B1=matrix(0,K,K)
  n.lk=matrix(0,K,K)
  N.lk=matrix(0,K,K)
  P0.avg=rep(0,max(stencil))
  llr.avg=0
  llr.rejection.avg=0
  t1.rejected.BH.chisq.avg=0
  t2.rejected.BH.chisq.avg=0
  llr.bic.avg=0
  
  left.right.param<-max(diag(stencil))
  
  anova.dfs=list()
  
  for(i in 1:max(stencil)){
    anova.dfs=append(anova.dfs, 
                     list(
                       'Data'=data.frame()
                     )
    )
  }
  
  comms=length(Tree.full$clct)
  
  left.right=list()
  for(i in 1:max(stencil)){
    left.right=append(left.right, 
                      list(
                        'Data'=list(mat=c())
                      )
    )
  }
  
  count=0
  
  for(i in 1:length(results.list.inner)){
    rejectedBH.avg=rejectedBH.avg+(results.list.inner[[i]]$rejectedBH.chisq)
    rejectedBF.avg=rejectedBF.avg+(results.list.inner[[i]]$rejectedBF.chisq)
    P0.avg=P0.avg+(results.list.inner[[i]]$P0)
    n.lk=n.lk+results.list.inner[[i]]$n.lk
    N.lk=N.lk+results.list.inner[[i]]$N.lk
    B0=B0+(results.list.inner[[i]]$B0)*results.list.inner[[i]]$n.lk
    B1=B1+(results.list.inner[[i]]$B1)*results.list.inner[[i]]$N.lk
    t1.rejected.BH.chisq.avg=t1.rejected.BH.chisq.avg+results.list.inner[[i]]$t1.rejectedBH.chisq
    t2.rejected.BH.chisq.avg=t2.rejected.BH.chisq.avg+results.list.inner[[i]]$t2.rejectedBH.chisq
    llr.rejection.avg=llr.rejection.avg+results.list.inner[[i]]$llr.test
    llr.bic.avg=llr.bic.avg+results.list.inner[[i]]$llr.bic
    for(j in 1:max(stencil)){
      B0.vec=results.list.inner[[i]]$B0[lower.tri(stencil,diag = T)]
      stencil.vec=stencil[lower.tri(stencil,diag = T)]
      vecP=B0.vec[stencil.vec==j]
      vecS=as.character(1:length(vecP))
      anova.dfs[[j]]=rbind(anova.dfs[[j]],data.frame('site'=vecS,'P'=vecP))
    }
    
    for(j in 1:max(stencil)){
      B0.vec=results.list.inner[[i]]$B0[lower.tri(stencil,diag = T)]
      stencil.vec=stencil[lower.tri(stencil,diag = T)]
      vec=B0.vec[stencil.vec==j]
      left.right[[j]]$mat=rbind(left.right[[j]]$mat,vec)
    }
    count=count+1
  }
  
  
  rejectedBH.avg=rejectedBH.avg/count
  rejectedBF.avg=rejectedBF.avg/count
  P0.avg=P0.avg/count
  BHline.avg=(1:length(P0.avg))*rep(0.05/length(P0.avg))
  B0=B0/n.lk
  B1=B1/N.lk
  t1.rejected.BH.chisq.avg=t1.rejected.BH.chisq.avg/count
  t2.rejected.BH.chisq.avg=t2.rejected.BH.chisq.avg/count
  llr.rejection.avg=llr.rejection.avg/count
  llr.bic.avg=llr.bic.avg/count
  
  
  P=unlist(lapply(anova.dfs, anova.test.lapply),use.names = F)
  P[which(is.nan(P))]=1
  anova.test=BH.rejection(P,stencil)
  
  rejectedBH.anova=anova.test$rejected
  
  t1.t2=t.error(rejectedBH.anova,pD,stencil)
  
  t1.rejected.BH.anova=t1.t2[1]
  t2.rejected.BH.anova=t1.t2[2]
  
  
  P=unlist(lapply(left.right, friedman.test.lapply.lr.mult),use.names = F)
  P[which(is.nan(P))]=1
  wilcox.lr=BH.rejection(P,stencil)
  
  rejected.wilcox.lr=wilcox.lr$rejected
  
  t1.t2=t.error(rejected.wilcox.lr,pD,stencil)
  
  t1.rejected.wilcox.lr=t1.t2[1]
  t2.rejected.wilcox.lr=t1.t2[2]
  
  gammaB=gamma.B(stencil)
  stat=2*LLRind(B0,B1,n.lk,stencil)
  P.chis.agg=pchisq(stat,gammaB,lower.tail = FALSE)
  P.chis.agg[is.nan(P.chis.agg)]=1
  P.chis.agg[which(B1==0)]=1
  P=c()
  
  for(i in 1:max(stencil)){
    P=c(P,P.chis.agg[which(stencil==i)][1])
  }
  
  BH.chisq.agg=BH.rejection(P,stencil)
  
  rejected.BH.chisq.agg=BH.chisq.agg$rejected
  
  t1.t2=t.error(rejected.BH.chisq.agg,pD,stencil)
  
  t1.rejected.BH.chisq.agg=t1.t2[1]
  t2.rejected.BH.chisq.agg=t1.t2[2]
  
  results.list.outer=append(results.list.outer,list(list(ind.results=results.list.inner,
                                                     org.stencil=org.stencil,
                                                     wilcox.lr=wilcox.lr,
                                                     anova.test=anova.test,
                                                     BH.chisq.agg=BH.chisq.agg,
                                                     rejectedBH.avg=rejectedBH.avg,
                                                     rejectedBF.avg=rejectedBF.avg,
                                                     P0.avg=P0.avg,
                                                     BHline=BHline.avg,
                                                     llr.rejection.avg=llr.rejection.avg,
                                                     llr.bic.avg=llr.bic.avg,
                                                     B0=B0,
                                                     B1=B1,
                                                     t1.rejected.BH.anova=t1.rejected.BH.anova,
                                                     t2.rejected.BH.anova=t2.rejected.BH.anova,
                                                     t1.rejected.wilcox.lr=t1.rejected.wilcox.lr,
                                                     t2.rejected.wilcox.lr=t2.rejected.wilcox.lr,
                                                     t1.rejected.BH.chisq.agg=t1.rejected.BH.chisq.agg,
                                                     t2.rejected.BH.chisq.agg=t2.rejected.BH.chisq.agg,
                                                     t1.rejected.BH.chisq.avg=t1.rejected.BH.chisq.avg,
                                                     t2.rejected.BH.chisq.avg=t2.rejected.BH.chisq.avg)))
  saveRDS(results.list.outer, file=FILESAVENAME)
}",sep=""), file=paste("runsim",i,".R",sep=""))}