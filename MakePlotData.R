sims=readRDS(file='models.RData')


big.df.BNU1.nocorr=data.frame()
big.df.BNU1.corr=data.frame()
big.df.3m.nocorr=data.frame()
big.df.3m.corr=data.frame()

plot.data=list("BNU1"=list(llr.test=data.frame(
                                               x=c(),
                                               y=c(),
                                               Model=c()
                                               ),
  
                          llr.BIC=data.frame(
                                               x=c(),
                                               y=c(),
                                               ymin=c(),
                                               ymax=c(),
                                               ymins=c(),
                                               ymaxs=c(),
                                               Model=c()
                                               ),
                          llr.BIC.lm=data.frame(x=c(),
                                                y=c(),
                                                Model=c()),
                          llr.BIC.multi=data.frame(x=c(),
                                                y=c(),
                                                Model=c(),
                                                run=c()),
                          resultscorr.chisq.avg=list(),
                          resultscorr.chisq.agg=list(),
                          resultscorr.friedman=list(),
                          resultscorr.anova=list(),
                          org.stencils=list()
                           ),
               "m3"= list(llr.test=data.frame(
                                              x=c(),
                                              y=c(),
                                              Model=c()
                                              ),
                          llr.BIC=data.frame(
                                              x=c(),
                                              y=c(),
                                              ymin=c(),
                                              ymax=c(),
                                              ymins=c(),
                                              ymaxs=c(),
                                              Model=c()
                                              ),
                          llr.BIC.lm=data.frame(x=c(),
                                                y=c(),
                                                Model=c()),
                          llr.BIC.multi=data.frame(x=c(),
                                                   y=c(),
                                                   Model=c(),
                                                   run=c()),
                          resultscorr.chisq.avg=list(),
                          resultscorr.chisq.agg=list(),
                          resultscorr.friedman=list(),
                          resultscorr.anova=list(),
                          org.stencils=list()
                          )
               )

nDs=c(0:10,seq(10,100,10))
for(corr in c("","-nocorr")){
  for(simnum in 1:21){
    for(model.name in c("model-brainnDs-","model-3m-7comnDs-")){#
      
      nD=nDs[simnum]
      
      file.name=paste(model.name,nD,corr,".Rdata",sep = "")
      
      if(file.exists(file.name)){
        tempdata=readRDS(file.name)  
        
        if(tempdata[[1]]$ind.results[[1]]$Tree=="3m-7com"){
          mod="m3"
        }else{
          mod="BNU1"
        }
        
        if(corr==""){
          M="With Variation"
        }else{
          M="No Variation"
        }
        
        org.stencil.temp=tempdata[[1]]$org.stencil
        
        plot.data[[mod]]$org.stencils[[paste(nD)]]=org.stencil.temp
        
        llr.BIC.temp=c()
        llr.test.temp=c()
        
        sten.dim=dim(org.stencil.temp)[1]
        
        P0.len=length(tempdata[[1]]$ind.results[[1]]$P0)
        
        count=sum(unlist(lapply(tempdata, function(x) return(length(x$ind.results)))))
        
        resultscorr.chisq.avg.temp=list("X"=data.frame(P0=rep(0,P0.len),n=rep(0,P0.len),r=rep(0,P0.len)),
                                        "rejected"=matrix(0,sten.dim,sten.dim))
        resultscorr.chisq.avg.temp.P0=c()
        
        resultscorr.chisq.agg.temp=list("X"=data.frame(P0=rep(0,P0.len),n=rep(0,P0.len),r=rep(0,P0.len)),
                                        "rejected"=matrix(0,sten.dim,sten.dim))
        resultscorr.chisq.agg.temp.P0=c()
        
        resultscorr.friedman.temp=list("X"=data.frame(P0=rep(0,P0.len),n=rep(0,P0.len),r=rep(0,P0.len)),
                                       "rejected"=matrix(0,sten.dim,sten.dim))
        resultscorr.friedman.temp.P0=c()
        
        resultscorr.anova.temp=list("X"=data.frame(P0=rep(0,P0.len),n=rep(0,P0.len),r=rep(0,P0.len)),
                                    "rejected"=matrix(0,sten.dim,sten.dim))
        resultscorr.anova.temp.P0=c()
        
        for(i in 1:length(tempdata)){
          llr.BIC.temp.temp=c()
          for(j in 1:length(tempdata[[i]]$ind.results)){
            llr.BIC.temp=c(llr.BIC.temp,tempdata[[i]]$ind.results[[j]]$llr.bic)
            llr.BIC.temp.temp=c(llr.BIC.temp.temp,tempdata[[i]]$ind.results[[j]]$llr.bic)
            llr.test.temp=c(llr.test.temp,tempdata[[i]]$ind.results[[j]]$llr.test)
            resultscorr.chisq.avg.temp$rejected=resultscorr.chisq.avg.temp$rejected+tempdata[[i]]$ind.results[[j]]$rejectedBH.chisq
            resultscorr.chisq.avg.temp$X=resultscorr.chisq.avg.temp$X+data.frame(P0=tempdata[[i]]$ind.results[[j]]$P0,
                                                                                 n=1:length(tempdata[[i]]$ind.results[[j]]$P0),
                                                                                 r=0.05*(1:length(tempdata[[i]]$ind.results[[j]]$P0))/length(tempdata[[i]]$ind.results[[j]]$P0)
                                                                                 )
            resultscorr.chisq.avg.temp.P0=rbind(resultscorr.chisq.avg.temp.P0,tempdata[[i]]$ind.results[[j]]$P0)
            
            
          }
          plot.data[[mod]]$llr.BIC.multi=rbind(plot.data[[mod]]$llr.BIC.multi, data.frame(
                                                                                          x=nD,
                                                                                          y=mean(llr.BIC.temp.temp),
                                                                                          Model=M,
                                                                                          run=paste("run",i,M,sep="")
                                                                                          )
                                               )
          
          resultscorr.chisq.agg.temp$rejected=resultscorr.chisq.agg.temp$rejected+tempdata[[i]]$BH.chisq.agg$rejected
          resultscorr.chisq.agg.temp$X=resultscorr.chisq.agg.temp$X+tempdata[[i]]$BH.chisq.agg$X
          resultscorr.chisq.agg.temp.P0=rbind(resultscorr.chisq.agg.temp.P0,tempdata[[i]]$BH.chisq.agg$X$P0)
          
          resultscorr.friedman.temp$rejected=resultscorr.friedman.temp$rejected+tempdata[[i]]$wilcox.lr$rejected
          resultscorr.friedman.temp$X=resultscorr.friedman.temp$X+tempdata[[i]]$wilcox.lr$X
          resultscorr.friedman.temp.P0=rbind(resultscorr.friedman.temp.P0,tempdata[[i]]$wilcox.lr$X$P0)
          
          resultscorr.anova.temp$rejected=resultscorr.anova.temp$rejected+tempdata[[i]]$anova.test$rejected
          resultscorr.anova.temp$X=resultscorr.anova.temp$X+tempdata[[i]]$anova.test$X
          resultscorr.anova.temp.P0=rbind(resultscorr.anova.temp.P0,tempdata[[i]]$anova.test$X$P0)
          print(paste(file.name,"iter",i,"step",j))
        }
        
        resultscorr.chisq.avg.temp$rejected=resultscorr.chisq.avg.temp$rejected/count
        resultscorr.chisq.avg.temp$X=resultscorr.chisq.avg.temp$X/count
        ymin=apply(resultscorr.chisq.avg.temp.P0,2,function(x) return(quantile(x,0.05 ,na.rm= TRUE)))
        ymax=apply(resultscorr.chisq.avg.temp.P0,2,function(x) return(quantile(x,0.95 ,na.rm= TRUE)))
        plot.data[[mod]]$resultscorr.chisq.avg[[paste(nD,corr)]]=list("X"=data.frame(P0=resultscorr.chisq.avg.temp$X$P0,
                                                                                n=resultscorr.chisq.avg.temp$X$n,
                                                                                r=resultscorr.chisq.avg.temp$X$r,
                                                                                P0min=ymin,
                                                                                P0max=ymax),
                                                                 "rejected"=resultscorr.chisq.avg.temp$rejected)
        
        resultscorr.chisq.agg.temp$rejected=resultscorr.chisq.agg.temp$rejected/length(tempdata)
        resultscorr.chisq.agg.temp$X=resultscorr.chisq.agg.temp$X/length(tempdata)
        ymin=apply(resultscorr.chisq.agg.temp.P0,2,function(x) return(quantile(x,0.05 ,na.rm= TRUE)))
        ymax=apply(resultscorr.chisq.agg.temp.P0,2,function(x) return(quantile(x,0.95 ,na.rm= TRUE)))
        plot.data[[mod]]$resultscorr.chisq.agg[[paste(nD,corr)]]=list("X"=data.frame(P0=resultscorr.chisq.agg.temp$X$P0,
                                                                                n=resultscorr.chisq.agg.temp$X$n,
                                                                                r=resultscorr.chisq.agg.temp$X$r,
                                                                                P0min=ymin,
                                                                                P0max=ymax),
                                                                 "rejected"=resultscorr.chisq.agg.temp$rejected)
      
        resultscorr.friedman.temp$rejected=resultscorr.friedman.temp$rejected/length(tempdata)
        resultscorr.friedman.temp$X=resultscorr.friedman.temp$X/length(tempdata)
        ymin=apply(resultscorr.friedman.temp.P0,2,function(x) return(quantile(x,0.05 ,na.rm= TRUE)))
        ymax=apply(resultscorr.friedman.temp.P0,2,function(x) return(quantile(x,0.95 ,na.rm= TRUE)))
        plot.data[[mod]]$resultscorr.friedman[[paste(nD,corr)]]=list("X"=data.frame(P0=resultscorr.friedman.temp$X$P0,
                                                                                n=resultscorr.friedman.temp$X$n,
                                                                                r=resultscorr.friedman.temp$X$r,
                                                                                P0min=ymin,
                                                                                P0max=ymax),
                                                                "rejected"=resultscorr.friedman.temp$rejected)
      
        resultscorr.anova.temp$rejected=resultscorr.anova.temp$rejected/length(tempdata)
        resultscorr.anova.temp$X=resultscorr.anova.temp$X/length(tempdata)
        ymin=apply(resultscorr.anova.temp.P0,2,function(x) return(quantile(x,0.05 ,na.rm= TRUE)))
        ymax=apply(resultscorr.anova.temp.P0,2,function(x) return(quantile(x,0.95 ,na.rm= TRUE)))
        plot.data[[mod]]$resultscorr.anova[[paste(nD,corr)]]=list("X"=data.frame(P0=resultscorr.anova.temp$X$P0,
                                                                            n=resultscorr.anova.temp$X$n,
                                                                            r=resultscorr.anova.temp$X$r,
                                                                            P0min=ymin,
                                                                            P0max=ymax),
                                                             "rejected"=resultscorr.anova.temp$rejected)
        
        
        
        
        
        plot.data[[mod]]$llr.BIC=rbind(plot.data[[mod]]$llr.BIC,data.frame(
                                                                              x=nD,
                                                                              y=mean(llr.BIC.temp ,na.rm= TRUE),
                                                                              ymin=quantile(llr.BIC.temp, probs = 0.05 ,na.rm= TRUE),
                                                                              ymax=quantile(llr.BIC.temp, probs = 0.95 ,na.rm= TRUE),
                                                                              ymins=mean(llr.BIC.temp ,na.rm= TRUE)-1.5*sqrt(var(llr.BIC.temp ,na.rm= TRUE)),
                                                                              ymaxs=mean(llr.BIC.temp ,na.rm= TRUE)+1.5*sqrt(var(llr.BIC.temp ,na.rm= TRUE)),
                                                                              Model=M
                                                                              )
                                        )
        
        plot.data[[mod]]$llr.BIC.lm=rbind(plot.data[[mod]]$llr.BIC.lm,data.frame(
                                                                                  x=nD,
                                                                                  y=llr.BIC.temp,
                                                                                  Model=M
                                                                                )
        )
        
        plot.data[[mod]]$llr.test=rbind(plot.data[[mod]]$llr.test,data.frame(
                                                                             x=nD,
                                                                             y=1-mean(llr.test.temp ,na.rm= TRUE),
                                                                             Model=M
                                                                             )
        )
      }
    }
  }
}

saveRDS(plot.data, "plotable_data_withlm_s.Rdata")








