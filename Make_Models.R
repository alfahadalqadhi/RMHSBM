

set.seed(1234)


nDs1=seq(0,10,1)
nDs2=seq(10,100,10)
nDs=c(nDs1,nDs2)




Tree.full.list=list("brain"=list(motif.sizes=c(35),motifs=c(1,1),clct=c(1,1),name="brain"),
                    "3m-7com"=list(motif.sizes=c(10,12,10),motifs=c(1,2,1,3,1,2,3),clct=c(2,1,2,2),name="3m-7com"))


models=list()

for(Tree.full in Tree.full.list){
  for(nD in nDs){
    results.list <- list()
    
    if(nD==0){
      org.stencil=get.stencil.general(Tree.full,nD=nD)
      pD=c()
    }else{
      org.stencil.pD=get.stencil.general(Tree.full,nD=nD)
      org.stencil=org.stencil.pD[[1]]
      pD=org.stencil.pD[[2]]
    }
    
    
    stencil=get.stencil.general(Tree.full,nD=0)
    
    models=append(models,list(list("nD"=nD,
                                   "model"=Tree.full,
                                   "org.stencil"=org.stencil,
                                   "pD"=pD)))
    
  }
}

saveRDS(models, file="models.RData")