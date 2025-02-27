
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



#' @export
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

#' @export
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


#' @export
LLRind <- function(B0hat, B1hat, n.lk, stencil) {
  K=dim(n.lk)[1]
  temp=n.lk * B0hat * log(B0hat / (1 - B0hat)) +
    n.lk * log(1 - B0hat) -
    n.lk * B0hat * log(B1hat / (1 - B1hat)) -
    n.lk * log(1 - B1hat)
  temp[which(is.nan(temp), arr.ind = T)]=0
  
  return(temp)
}

#' @export
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

#' @export

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

#' @export

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

#' @export

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

#' @export
#' 

make.graph = function(n, B, block.sizes) {
  graph_sbm = sample_sbm(n,
                         B,
                         block.sizes = block.sizes,
                         directed = FALSE,
                         loops = FALSE)
  A = matrix(as_adjacency_matrix(graph_sbm), nrow = n)
  return(A)
}

#' @export

get.t = function(vec) {
  temp = c()
  for (i in 1:length(vec)) {
    temp = c(temp, rep(i, vec[i]))
  }
  return(temp)
}


#' @export

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

#' @export

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
  return(summary(aov(P~site,x))[[1]][["Pr(>F)"]][1])
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
  
  return(list("rejected"=rejected, "X"=X))
}

friedman.test.lapply.lr.mult= function(x){
  return(friedman.test(x$mat)$p.value)
}

plot.GG.blocks.general.color.stencil= function(X,title,Tree.full,org.stencil){
  stencil=get.stencil.general(Tree.full,nD=0)
  comm.sizes=Tree.full$motif.sizes[Tree.full$motifs]
  G=1:sum(comm.sizes)
  clct=Tree.full$clct
  mx=max(stencil)
  X[which(org.stencil>mx)]=-1
  X.df <-
    X[G,G]%>%
    as_tibble() %>%
    add_column("Var1"=G,.before = 1) %>%
    rename_all(~tibble(col=c("Var1",paste("V",G,sep = ""))) %>% pull(col)) %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
    mutate(
      Var1 = factor(Var1, levels = G),
      Var2 = factor(gsub("V", "", Var2), levels = rev(G))
    )%>%
    add_column("Zero"=as.integer(as.logical(org.stencil>mx)),.after = 3)
  
  p0=ggplot(X.df, aes(Var1, Var2)) +
    geom_tile(aes(fill = value))+
    scale_fill_gradientn(colours = c("#295F98","white","black"),
                         values = rescale(c(-1,0,mx)),
                         limits=c(-1,mx))+
    geom_segment(data=data.frame(x=0.5,
                                 y=0.5,
                                 xend=0.5,
                                 yend=max(G)+0.5),
                 aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
    geom_segment(data=data.frame(x=0.5,
                                 y=max(G)+0.5,
                                 xend=max(G)+0.5,
                                 yend=max(G)+0.5),
                 aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
    xlab("")+ylab("")+ggtitle(title)+
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust=0.5),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  passed=0
  count=1
  passed.clct=0
  for(i in comm.sizes){
    passed=passed+i
    passed.clct=passed.clct+i
    clct.end=sum(comm.sizes[1:sum(clct[1:count])])
    if(passed<clct.end){
      p0<-p0+geom_segment(data=data.frame(x=passed+0.5,
                                          y=max(G)-clct.end+0.5,
                                          xend=passed+0.5,
                                          yend=max(G)-passed+passed.clct+0.5),
                          aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
        geom_segment(data=data.frame(x=passed-passed.clct+0.5,
                                     y=max(G)-passed+0.5,
                                     xend=clct.end+0.5,
                                     yend=max(G)-passed+0.5),
                     aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")
      
    }else{
      p0=p0+geom_segment(data=data.frame(x=passed+0.5,
                                         y=0.5,
                                         xend=passed+0.5,
                                         yend=max(G)+0.5),
                         aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
        geom_segment(data=data.frame(x=0.5,
                                     y=max(G)-passed+0.5,
                                     xend=max(G)+0.5,
                                     yend=max(G)-passed+0.5),
                     aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")
      count=count+1
      passed.clct=0
    }
    
    
    
  }
  return(p0)
}


plot.GG.blocks.general.color10= function(X,title,Tree.full,org.stencil){
  stencil=get.stencil.general(Tree.full,nD=0)
  comm.sizes=Tree.full$motif.sizes[Tree.full$motifs]
  G=1:sum(comm.sizes)
  clct=Tree.full$clct
  mx=max(stencil)
  X[which(org.stencil<=mx)]=-10*X[which(org.stencil<=mx)]
  X[which(org.stencil>mx)]=X[which(org.stencil>mx)]+1
  X.df <-
    X[G,G]%>%
    as_tibble() %>%
    add_column("Var1"=G,.before = 1) %>%
    rename_all(~tibble(col=c("Var1",paste("V",G,sep = ""))) %>% pull(col)) %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
    mutate(
      Var1 = factor(Var1, levels = G),
      Var2 = factor(gsub("V", "", Var2), levels = rev(G))
    )%>%
    add_column("Zero"=as.integer(as.logical(org.stencil>mx)),.after = 3)
  
  p0=ggplot(X.df, aes(Var1, Var2)) +
    geom_tile(aes(fill = value))+
    scale_fill_gradientn(colours = c("#0b0b60","white","#C96868","black"),
                         values = rescale(c(-1,0,1,2)),
                         limits=c(-1,2))+
    geom_segment(data=data.frame(x=0.5,
                                 y=0.5,
                                 xend=0.5,
                                 yend=max(G)+0.5),
                 aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
    geom_segment(data=data.frame(x=0.5,
                                 y=max(G)+0.5,
                                 xend=max(G)+0.5,
                                 yend=max(G)+0.5),
                 aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
    xlab("")+ylab("")+ggtitle(title)+
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust=0.5),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  passed=0
  count=1
  passed.clct=0
  for(i in comm.sizes){
    passed=passed+i
    passed.clct=passed.clct+i
    clct.end=sum(comm.sizes[1:sum(clct[1:count])])
    if(passed<clct.end){
      p0<-p0+geom_segment(data=data.frame(x=passed+0.5,
                                          y=max(G)-clct.end+0.5,
                                          xend=passed+0.5,
                                          yend=max(G)-passed+passed.clct+0.5),
                          aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
        geom_segment(data=data.frame(x=passed-passed.clct+0.5,
                                     y=max(G)-passed+0.5,
                                     xend=clct.end+0.5,
                                     yend=max(G)-passed+0.5),
                     aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")
      
    }else{
      p0=p0+geom_segment(data=data.frame(x=passed+0.5,
                                         y=0.5,
                                         xend=passed+0.5,
                                         yend=max(G)+0.5),
                         aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
        geom_segment(data=data.frame(x=0.5,
                                     y=max(G)-passed+0.5,
                                     xend=max(G)+0.5,
                                     yend=max(G)-passed+0.5),
                     aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")
      count=count+1
      passed.clct=0
    }
    
    
    
  }
  return(p0)
}

plot.GG.blocks.general.color= function(X,title,Tree.full,org.stencil){
  stencil=get.stencil.general(Tree.full,nD=0)
  comm.sizes=Tree.full$motif.sizes[Tree.full$motifs]
  G=1:sum(comm.sizes)
  clct=Tree.full$clct
  mx=max(stencil)
  X[which(org.stencil<=mx)]=-X[which(org.stencil<=mx)]
  X[which(org.stencil>mx)]=X[which(org.stencil>mx)]+1
  X.df <-
    X[G,G]%>%
    as_tibble() %>%
    add_column("Var1"=G,.before = 1) %>%
    rename_all(~tibble(col=c("Var1",paste("V",G,sep = ""))) %>% pull(col)) %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
    mutate(
      Var1 = factor(Var1, levels = G),
      Var2 = factor(gsub("V", "", Var2), levels = rev(G))
    )%>%
    add_column("Zero"=as.integer(as.logical(org.stencil>mx)),.after = 3)
  
  p0=ggplot(X.df, aes(Var1, Var2)) +
    geom_tile(aes(fill = value))+
    scale_fill_gradientn(colours = c("#0b0b60","white","#C96868","black"),
                         values = rescale(c(-1,0,1,2)),
                         limits=c(-1,2))+
    geom_segment(data=data.frame(x=0.5,
                                 y=0.5,
                                 xend=0.5,
                                 yend=max(G)+0.5),
                 aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
    geom_segment(data=data.frame(x=0.5,
                                 y=max(G)+0.5,
                                 xend=max(G)+0.5,
                                 yend=max(G)+0.5),
                 aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
    xlab("")+ylab("")+ggtitle(title)+
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust=0.5),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  passed=0
  count=1
  passed.clct=0
  for(i in comm.sizes){
    passed=passed+i
    passed.clct=passed.clct+i
    clct.end=sum(comm.sizes[1:sum(clct[1:count])])
    if(passed<clct.end){
      p0<-p0+geom_segment(data=data.frame(x=passed+0.5,
                                          y=max(G)-clct.end+0.5,
                                          xend=passed+0.5,
                                          yend=max(G)-passed+passed.clct+0.5),
                          aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
        geom_segment(data=data.frame(x=passed-passed.clct+0.5,
                                     y=max(G)-passed+0.5,
                                     xend=clct.end+0.5,
                                     yend=max(G)-passed+0.5),
                     aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")
      
    }else{
      p0=p0+geom_segment(data=data.frame(x=passed+0.5,
                                         y=0.5,
                                         xend=passed+0.5,
                                         yend=max(G)+0.5),
                         aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")+
        geom_segment(data=data.frame(x=0.5,
                                     y=max(G)-passed+0.5,
                                     xend=max(G)+0.5,
                                     yend=max(G)-passed+0.5),
                     aes(x=x,y=y,xend=xend,yend=yend),colour="grey30")
      count=count+1
      passed.clct=0
    }
    
    
    
  }
  return(p0)
}

plot.P=function(X){
  return(ggplot(X,aes(x=n,y=P0,ymin=P0min,ymax=P0max))+
         geom_line()+geom_ribbon(alpha=0.5)+geom_line(aes(x=n,y=r), linetype="dashed",lwd=0.8)+
         ggtitle("")+theme(plot.title = element_text(hjust=0.5))+
         xlab("Order")+ylab("P"))
}

plot.P.noribbon=function(X){
  return(ggplot(X,aes(x=n,y=P0))+
           geom_line()+geom_line(aes(x=n,y=r), linetype="dashed",lwd=0.8)+
           ggtitle("")+theme(plot.title = element_text(hjust=0.5))+
           xlab("Order")+ylab("P"))
}
