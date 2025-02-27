BNU1=readRDS(file="DataBNU1_RMHSBM6_BH_BF_P.RData")
subjects=1:57
scan=c(1)
ni=c(300)

rejectedBH.avg=matrix(0,70,70)
rejectedBF.avg=matrix(0,70,70)
B0=matrix(0,70,70)
B1=matrix(0,70,70)
n.lk=matrix(0,70,70)
N.lk=matrix(0,70,70)
P0.avg=rep(0,631)
BHline.avg=rep(0,631)

stencil=get.stencil(70)

left.right=list()
for(i in 1:max(stencil)){
  left.right=append(left.right, 
                    list(
                      "Data"=list(mat=c())
                    )
  )
}



anova.dfs=list()

for(i in 1:max(stencil)){
  anova.dfs=append(anova.dfs, 
                   list(
                     "Data"=data.frame()
                   )
  )
}

count=0

for(i in 1:length(BNU1)){
  if(BNU1[[i]]$Subject %in% subjects){
    if(BNU1[[i]]$Scan %in% scan){
      if(BNU1[[i]]$'n*' %in% ni){
        rejectedBH.avg=rejectedBH.avg+(BNU1[[i]]$rejectedBH)
        rejectedBF.avg=rejectedBF.avg+(BNU1[[i]]$rejectedBF)
        P0.avg=P0.avg+(BNU1[[i]]$P0)
        BHline.avg=BHline.avg+(BNU1[[i]]$BHline)
        n.lk=n.lk+BNU1[[i]]$n.lk
        N.lk=N.lk+BNU1[[i]]$N.lk
        B0=B0+(BNU1[[i]]$B0)*BNU1[[i]]$n.lk
        B1=B1+(BNU1[[i]]$B1)*BNU1[[i]]$N.lk
        
        
        for(j in 1:max(stencil)){
          B0.vec=BNU1[[i]]$B0[lower.tri(stencil,diag = T)]
          stencil.vec=stencil[lower.tri(stencil,diag = T)]
          vecP=B0.vec[stencil.vec==j]
          vecS=as.character(1:length(vecP))
          anova.dfs[[j]]=rbind(anova.dfs[[j]],data.frame("site"=vecS,"P"=vecP))
        }
        
        for(j in 1:max(stencil)){
          B0.vec=BNU1[[i]]$B0[lower.tri(stencil,diag = T)]
          stencil.vec=stencil[lower.tri(stencil,diag = T)]
          vec=B0.vec[stencil.vec==j]
          left.right[[j]]$mat=rbind(left.right[[j]]$mat,vec)
        }
        
        count=count+1
      }
    }
  }
}



rejectedBH.avg=rejectedBH.avg/count
rejectedBF.avg=rejectedBF.avg/count
P0.avg=P0.avg/count
BHline.avg=BHline.avg/count
B0=B0/n.lk
B1=B1/N.lk

P=unlist(lapply(anova.dfs, anova.test.lapply),use.names = F)
P[which(is.nan(P))]=1
anova.test=BH.rejection(P,stencil)

P=unlist(lapply(left.right, friedman.test.lapply.lr.mult),use.names = F)
P[which(is.nan(P))]=1
wilcox.lr=BH.rejection(P,stencil)

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

G=2:71

sdf=1:max(stencil)

X.chisq=data.frame(P0=P0.avg,n=sdf,r=BHline.avg)


plot.GG.blocks.general.cross= function(X,title,Tree.full,org.stencil){
  stencil=get.stencil.general(Tree.full,nD=0)
  comm.sizes=Tree.full$motif.sizes[Tree.full$motifs]
  G=1:sum(comm.sizes)
  clct=Tree.full$clct
  mx=max(stencil)
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
    scale_fill_gradient(low = "white", high = "black")+
    geom_point(data=X.df[which(X.df$Zero>0),],
               aes(Var1,Var2),
               shape=21,colour="black",
               fill="white",
               stroke=1.1,
               size=1.25)+
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
    scale_fill_gradientn(colours = c("black","white","#C96868","black"),
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



plot.result=function(BH.chisq.ind,BH.chisq.agg,anova.test,wilcox.lr,G,title,B1){
  stencil.8=get.stencil.general(Tree.full,nD=0)
  mx=max(stencil)
  stencil.8[which(B1==0)]=mx+1
  tree=list(motif.sizes=c(35),motifs=c(1,1),clct=c(1,1),name="brain")
  
  p0=plot.GG.blocks.general.color(BH.chisq.ind$rejected,"Chi−square averaged over data",tree,stencil.8)
  p1=plot.GG.blocks.general.color(BH.chisq.agg$rejected,"Chi−square on aggregated data",tree,stencil.8)
  p2=plot.GG.blocks.general.color(anova.test$rejected,"ANOVA",tree,stencil.8)
  p3=plot.GG.blocks.general.color(wilcox.lr$rejected,"Friedman",tree,stencil.8)
  
  
  p4=ggplot(BH.chisq.ind$X,aes(x=n,y=P0))+
    geom_line()+geom_line(aes(x=n,y=r), linetype="dashed",lwd=0.8)+
    ggtitle("P profile")+theme(plot.title = element_text(hjust=0.5))+
    xlab("Order")+ylab("P")
  p5=ggplot(BH.chisq.agg$X,aes(x=n,y=P0))+
    geom_line()+geom_line(aes(x=n,y=r), linetype="dashed",lwd=0.8)+
    ggtitle("P profile")+theme(plot.title = element_text(hjust=0.5))+
    xlab("Order")+ylab("P")
  p6=ggplot(anova.test$X,aes(x=n,y=P0))+
    geom_line()+geom_line(aes(x=n,y=r), linetype="dashed",lwd=0.8)+
    ggtitle("P profile")+theme(plot.title = element_text(hjust=0.5))+
    xlab("Order")+ylab("P")
  p7=ggplot(wilcox.lr$X,aes(x=n,y=P0))+
    geom_line()+geom_line(aes(x=n,y=r), linetype="dashed",lwd=0.8)+
    ggtitle("P profile")+theme(plot.title = element_text(hjust=0.5))+
    xlab("Order")+ylab("P")
  return(grid.arrange(p0,p1,p2,p3,p4,p5,p6,p7, nrow = 2))
}

BH.chisq.ind=list(X=X.chisq, rejected=rejectedBH.avg)

Tree.full=list(motif.sizes=c(35),motifs=c(1,1),clct=c(1,1),name="brain")
pdf(file="BNU1dataresults.pdf",
    width = 28,
    height = 14)
plot.result(BH.chisq.ind,BH.chisq.agg,anova.test,wilcox.lr,G,"ANOVA",B1)
dev.off()
