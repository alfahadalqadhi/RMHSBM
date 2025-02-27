

n=300



results.list <- list()

for(subj in c(1:57)){
  for(scan in c(1,2)){
    for(ni in n){
      # file=sprintf("https://www.cis.jhu.edu/~parky/Microsoft/JHU-MSR/ZMx2/BNU1/DS16784/subj%d-scan%d.graphml",
      #              subj,scan)
      
      file=sprintf("C:/Users/alfah/Documents/BNU116K/subj%d-scan%d.graphml",
                   subj,scan)
      
      
      network<-read_graph(file=file,
                          format = "graphml")
      
      block_memberships<-vertex_attr(network,name="region",index=V(network))
      
      G=(1:71)
      
      A=as_adjacency_matrix(network,sparse = F)
      
      B0=matrix(0,length(G),length(G))
      n.lk=matrix(0,length(G),length(G))
      
      for(l in 1:length(G)){
        for(k in 1:l){
          temp=A[block_memberships==G[l],block_memberships==G[k]]
          edges=sum(temp)
          size=length(c(temp))
          if(l==k){
            n.lk[l,k]= choose(sum(block_memberships==G[l]),2)
            B0[l,k]= edges/(2*n.lk[l,k])
          }else{
            n.lk[l,k]=n.lk[k,l]=size
            B0[l,k]=B0[k,l]=edges/size
          }
        }
      }
      
      too.smallh=c(unique(c(which(n.lk<= ni, arr.ind = T))))
      too.small1.pairs=too.smallh[which(too.smallh<=36)]+35
      too.small2.pairs=too.smallh[which(too.smallh>36)]-35
      too.small=unique(c(too.smallh,too.small1.pairs,too.small2.pairs))
      
      too.different=unique(c(which(abs(B0[2:36,2:36]-B0[37:71,37:71])>=0.018,arr.ind = T)+1))
      
      #exclusions=unique(c(too.small,too.different,too.different+35,1))
      exclusions<-c(1)
      G=(1:71)[! (1:71) %in% exclusions]
      # G=c(c(2,3,9,10,11,12,25,26,27,31,32,33),c(2,3,9,10,11,12,25,26,27,31,32,33)+35)
      B0=matrix(0,length(G),length(G))
      n.lk=matrix(0,length(G),length(G))
      
      for(l in 1:length(G)){
        for(k in 1:l){
          temp=A[block_memberships==G[l],block_memberships==G[k]]
          edges=sum(temp)
          size=length(c(temp))
          if(l==k){
            n.lk[l,k]= choose(sum(block_memberships==G[l]),2)
            B0[l,k]= edges/(2*n.lk[l,k])
          }else{
            n.lk[l,k]=n.lk[k,l]=size
            B0[l,k]=B0[k,l]=edges/size
          }
        }
      }
      
      
      comm.size=length(G)/2
      
      stencil_half=matrix(0,comm.size,comm.size)
      stencil_half[lower.tri(stencil_half,diag = T)]=1:(choose(comm.size,2)+comm.size)
      lowerstencil=matrix((choose(comm.size,2)+comm.size+1),comm.size,comm.size)
      upperstencil=matrix(0,comm.size,comm.size)
      stenciltop=cbind(stencil_half,upperstencil)
      stencilbottom=cbind(lowerstencil,stencil_half)
      stencil=rbind(stenciltop,stencilbottom)
      
      N.lk=matrix(0,length(G),length(G))
      B1.diag=rep(0,length(G))
      
      
      for (l in unique(diag(stencil))){
        blocks=which(stencil==l,arr.ind = T)[,1]
        edges=0
        size=0
        for (b in blocks){
          temp=A[block_memberships==G[b],block_memberships==G[b]]
          edges=edges+sum(temp)/2
          size=size+n.lk[b,b]
        }
        B1.diag[blocks]=edges/size
        diag(N.lk)[blocks]=size
      }
      
      B1=matrix(0,length(G),length(G))
      for (l in unique(stencil[lower.tri(stencil)])){
        blocks.1=which(stencil==l,arr.ind = T)[,2]
        blocks.2=which(stencil==l,arr.ind = T)[,1]
        edges=0
        size=0
        for (b in 1:length(blocks.1)){
          temp=A[block_memberships==G[blocks.1[b]],block_memberships==G[blocks.2[b]]]
          edges=edges+sum(temp)
          size=size+length(c(temp))
          B1[blocks.1[b],blocks.2[b]]=B1[blocks.2[b],blocks.1[b]]=NA
        }
        N.lk[which(is.na(B1))]=size
        B1[which(is.na(B1))]=edges/size
      }
      
      diag(B1)=B1.diag
      


      gammaB=gamma.B(stencil)
      stat=2*LLRind(B0,B1,n.lk,stencil)
      P=pchisq(stat,gammaB,lower.tail = FALSE)
      P[upper.tri(P)]=t(P)[upper.tri(P)]
      P[is.nan(P)]=1
      P[which(B0==0)]=1
      P0=c()
      
      for(i in 1:max(stencil)){
        P0=c(P0,P[which(stencil==i)][1])
      }
      P0.ordered.chisq=P0[order(P0)]
      param=1:max(stencil)
      param.ordered=param[order(P0)]
      regline=(1:length(P0))*rep(0.05/length(P0))
      k=max(which(P0.ordered.chisq<regline))
      if(k!=0){
        ind=which(stencil %in% param.ordered[1:k])
      }else{
        ind=c()
      }
      
      
      rejectedBH.chisq=matrix(0,length(G),length(G))
      rejectedBH.chisq[ind]=1
      rejectedBH.chisq[upper.tri(P)]=t(rejectedBH.chisq)[upper.tri(P)]
      
      
      bonferroni=qchisq(0.05/(max(stencil)),gammaB)
      rejectedBF=matrix(0,length(G),length(G))
      rejectedBF[which(P<bonferroni)]=1
      rejectedBF[upper.tri(P)]=t(rejectedBF)[upper.tri(P)]
      
      sdf=1:length(P0)
      
      range=0:100
      
      X=data.frame(n=1:length(P0), P0=P0.ordered.chisq,r=regline)
      
           
      grid.arrange(p0,p1,p2,p3, nrow = 2)
        
      
      results.list <- append(results.list, list("Data"=list(
        "Subject"=subj,
        "Scan"=scan,
        "n*"=ni,
        "B0"=B0,
        "B1"=B1,
        "n.lk"=n.lk,
        "N.lk"=N.lk,
        "rejectedBH.chisq"=rejectedBH.chisq,
        "rejectedBF"=rejectedBF,
        "P0"=P0.ordered.chisq,
        "BHline"=regline
      )))
      
      print(sprintf("subj: %f, scan: %f, n*=%d", subj, scan,ni))
      
    }
  }
}


saveRDS(results.list, file="DataBNU1_RMHSBM6_BH_BF_P.RData")