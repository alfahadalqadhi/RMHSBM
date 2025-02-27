
BNU1=70
m3=70

plot.data=readRDS("plotable_data_withlm_s.Rdata")

Tree.full.list=list("brain"=list(motif.sizes=c(35),motifs=c(1,1),clct=c(1,1),name="brain"),
                    "3m-7com"=list(motif.sizes=c(10,12,10),motifs=c(1,2,1,3,1,2,3),clct=c(2,1,2,2),name="3m-7com"))

stencil.plot.BNU1=get.stencil.general(Tree.full.list[["brain"]])

org.stencil.BNU1=plot.data[["BNU1"]]$org.stencils[[paste(BNU1)]]

mx=max(stencil.plot.BNU1)
mxo=max(org.stencil.BNU1)

for(i in mx:mxo){
  stencil.plot.BNU1[which(stencil.plot.BNU1==stencil.plot.BNU1[which(org.stencil.BNU1==i)][1])]=i
}


stencil.plot.m3=get.stencil.general(Tree.full.list[["3m-7com"]])



org.stencil.m3=plot.data[["m3"]]$org.stencils[[paste(m3)]]

mx=max(stencil.plot.m3)
mxo=max(org.stencil.m3)

for(i in mx:mxo){
  stencil.plot.m3[which(stencil.plot.m3==stencil.plot.m3[which(org.stencil.m3==i)][1])]=i
}

pBNU1stencil1=plot.GG.blocks.general.color.stencil(org.stencil.BNU1,"",Tree.full.list[["brain"]],stencil.plot.BNU1)
p3mstencil1=plot.GG.blocks.general.color.stencil(org.stencil.m3,"",Tree.full.list[["3m-7com"]],stencil.plot.m3)

pdf(file="BNU1results.pdf",
    width = 28,
    height = 14)


p0=plot.GG.blocks.general.color10(plot.data[["BNU1"]]$resultscorr.chisq.avg[[paste(BNU1,"-nocorr")]]$rejected,
                                  "Chi-square averaged over data",
                                  Tree.full.list[["brain"]],
                                  stencil.plot.BNU1)

p1=plot.GG.blocks.general.color10(plot.data[["BNU1"]]$resultscorr.chisq.agg[[paste(BNU1,"-nocorr")]]$rejected,
                                  "Chi-square on aggregated data",
                                  Tree.full.list[["brain"]],
                                  stencil.plot.BNU1)

p2=plot.GG.blocks.general.color10(plot.data[["BNU1"]]$resultscorr.anova[[paste(BNU1,"-nocorr")]]$rejected,
                                  "ANOVA",
                                  Tree.full.list[["brain"]],
                                  stencil.plot.BNU1)

p3=plot.GG.blocks.general.color10(plot.data[["BNU1"]]$resultscorr.friedman[[paste(BNU1,"-nocorr")]]$rejected,
                                  "Friedman",
                                  Tree.full.list[["brain"]],
                                  stencil.plot.BNU1)


p4=plot.P(plot.data[["BNU1"]]$resultscorr.chisq.avg[[paste(BNU1,"-nocorr")]]$X)
p5=plot.P(plot.data[["BNU1"]]$resultscorr.chisq.agg[[paste(BNU1,"-nocorr")]]$X)
p6=plot.P(plot.data[["BNU1"]]$resultscorr.anova[[paste(BNU1,"-nocorr")]]$X)
p7=plot.P(plot.data[["BNU1"]]$resultscorr.friedman[[paste(BNU1,"-nocorr")]]$X)

grid.arrange(p0,p1,p2,p3,p4,p5,p6,p7, nrow = 2)

dev.off()

pdf(file="BNU1resultscorr.pdf",
    width = 28,
    height = 14)


p0=plot.GG.blocks.general.color(plot.data[["BNU1"]]$resultscorr.chisq.avg[[paste(BNU1, "")]]$rejected,
                                "Chi-square averaged over data",
                                Tree.full.list[["brain"]],
                                stencil.plot.BNU1)

p1=plot.GG.blocks.general.color(plot.data[["BNU1"]]$resultscorr.chisq.agg[[paste(BNU1, "")]]$rejected,
                                "Chi-square on aggregated data",
                                Tree.full.list[["brain"]],
                                stencil.plot.BNU1)

p2=plot.GG.blocks.general.color10(plot.data[["BNU1"]]$resultscorr.anova[[paste(BNU1, "")]]$rejected,
                                "ANOVA",
                                Tree.full.list[["brain"]],
                                stencil.plot.BNU1)

p3=plot.GG.blocks.general.color10(plot.data[["BNU1"]]$resultscorr.friedman[[paste(BNU1, "")]]$rejected,
                                "Friedman",
                                Tree.full.list[["brain"]],
                                stencil.plot.BNU1)


p4=plot.P(plot.data[["BNU1"]]$resultscorr.chisq.avg[[paste(BNU1, "")]]$X)
p5=plot.P(plot.data[["BNU1"]]$resultscorr.chisq.agg[[paste(BNU1, "")]]$X)
p6=plot.P(plot.data[["BNU1"]]$resultscorr.anova[[paste(BNU1, "")]]$X)
p7=plot.P(plot.data[["BNU1"]]$resultscorr.friedman[[paste(BNU1, "")]]$X)

grid.arrange(p0,p1,p2,p3,p4,p5,p6,p7, nrow = 2)

dev.off()

pdf(file="3mresults.pdf",
    width = 28,
    height = 14)

p0=plot.GG.blocks.general.color10(plot.data[["m3"]]$resultscorr.chisq.avg[[paste(m3,"-nocorr")]]$rejected,
                                  "Chi-square averaged over data",
                                  Tree.full.list[["3m-7com"]],
                                  stencil.plot.m3)

p1=plot.GG.blocks.general.color10(plot.data[["m3"]]$resultscorr.chisq.agg[[paste(m3,"-nocorr")]]$rejected,
                                  "Chi-square on aggregated data",
                                  Tree.full.list[["3m-7com"]],
                                  stencil.plot.m3)

p2=plot.GG.blocks.general.color10(plot.data[["m3"]]$resultscorr.anova[[paste(m3,"-nocorr")]]$rejected,
                                  "ANOVA",
                                  Tree.full.list[["3m-7com"]],
                                  stencil.plot.m3)

p3=plot.GG.blocks.general.color10(plot.data[["m3"]]$resultscorr.friedman[[paste(m3,"-nocorr")]]$rejected,
                                  "Friedman",
                                  Tree.full.list[["3m-7com"]],
                                  stencil.plot.m3)


p4=plot.P(plot.data[["m3"]]$resultscorr.chisq.avg[[paste(m3,"-nocorr")]]$X)
p5=plot.P(plot.data[["m3"]]$resultscorr.chisq.agg[[paste(m3,"-nocorr")]]$X)
p6=plot.P(plot.data[["m3"]]$resultscorr.anova[[paste(m3,"-nocorr")]]$X)
p7=plot.P(plot.data[["m3"]]$resultscorr.friedman[[paste(m3,"-nocorr")]]$X)

grid.arrange(p0,p1,p2,p3,p4,p5,p6,p7, nrow = 2)

dev.off()


pdf(file="3mresultscorr.pdf",
    width = 28,
    height = 14)

p0=plot.GG.blocks.general.color(plot.data[["m3"]]$resultscorr.chisq.avg[[paste(m3,"")]]$rejected,
                                "Chi-square averaged over data",
                                Tree.full.list[["3m-7com"]],
                                stencil.plot.m3)

p1=plot.GG.blocks.general.color(plot.data[["m3"]]$resultscorr.chisq.agg[[paste(m3,"")]]$rejected,
                                "Chi-square on aggregated data",
                                Tree.full.list[["3m-7com"]],
                                stencil.plot.m3)

p2=plot.GG.blocks.general.color10(plot.data[["m3"]]$resultscorr.anova[[paste(m3,"")]]$rejected,
                                "ANOVA",
                                Tree.full.list[["3m-7com"]],
                                stencil.plot.m3)

p3=plot.GG.blocks.general.color10(plot.data[["m3"]]$resultscorr.friedman[[paste(m3,"")]]$rejected,
                                "Friedman",
                                Tree.full.list[["3m-7com"]],
                                stencil.plot.m3)


p4=plot.P(plot.data[["m3"]]$resultscorr.chisq.avg[[paste(m3,"")]]$X)
p5=plot.P(plot.data[["m3"]]$resultscorr.chisq.agg[[paste(m3,"")]]$X)
p6=plot.P(plot.data[["m3"]]$resultscorr.anova[[paste(m3,"")]]$X)
p7=plot.P(plot.data[["m3"]]$resultscorr.friedman[[paste(m3,"")]]$X)

grid.arrange(p0,p1,p2,p3,p4,p5,p6,p7, nrow = 2)

dev.off()




pdf(file="simsLLR.pdf",
    width = 16,
    height = 7)


scaleFactor0 <- max(abs(plot.data[["BNU1"]]$llr.BIC$ymax)) / max(plot.data[["BNU1"]]$llr.test$y)
p0=ggplot(plot.data[["BNU1"]]$llr.BIC, aes(x=x, y=y, ymin=ymin, ymax=ymax, fill=Model)) + 
  geom_line(lwd=1) + 
  geom_ribbon(alpha=0.5)+ggtitle("BNU1 Simulation")+
  xlab("# of Mismatching Parameters")+
  geom_line(aes(x=plot.data[["BNU1"]]$llr.test$x,
                y=plot.data[["BNU1"]]$llr.test$y*scaleFactor0,
                color=plot.data[["BNU1"]]$llr.test$Model),lwd=1, linetype=2)+
  geom_hline(yintercept = 0,lwd=1)+
  scale_y_continuous(name="BIC Penalized LLR", sec.axis=sec_axis(~./scaleFactor0, name="Rejection %")) +
  theme(
    plot.title = element_text(hjust=0.5),
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"),
    axis.line.y.left = element_line(color = "black",
                             linewidth = 1,
                             linetype = 1),
    axis.title.y.right=element_text(color="black"),
    axis.text.y.right=element_text(color="black"),
    axis.line.y.right = element_line(color = "black",
                                    linewidth = 1,
                                    linetype = 2),
    legend.position = c(0.85,0.1)
  )+labs(linetype="", colour="",fill="")


scaleFactor1 <- max(abs(plot.data[["m3"]]$llr.BIC$ymax)) / max(plot.data[["BNU1"]]$llr.test$y)
p1=ggplot(plot.data[["m3"]]$llr.BIC, aes(x=x, y=y, ymin=ymin, ymax=ymax, fill=Model)) + 
  geom_line(lwd=1) + 
  geom_ribbon(alpha=0.5)+ggtitle("3 Motif 7 Communities Simulation")+
  xlab("# of Mismatching Parameters")+
  geom_line(aes(x=plot.data[["m3"]]$llr.test$x,
                y=plot.data[["m3"]]$llr.test$y*scaleFactor1,
                color=plot.data[["m3"]]$llr.test$Model),lwd=1, linetype=2)+
  geom_hline(yintercept = 0,lwd=1)+
  scale_y_continuous(name="BIC Penalized LLR", sec.axis=sec_axis(~./scaleFactor0, name="Rejection %")) +
  theme(
    plot.title = element_text(hjust=0.5),
    axis.title.y.left=element_text(color="black"),
    axis.text.y.left=element_text(color="black"),
    axis.line.y.left = element_line(color = "black",
                                    linewidth = 1,
                                    linetype = 1),
    axis.title.y.right=element_text(color="black"),
    axis.text.y.right=element_text(color="black"),
    axis.line.y.right = element_line(color = "black",
                                     linewidth = 1,
                                     linetype = 2),
    legend.position = c(0.85,0.1)
  )+labs(linetype="", colour="",fill="")
grid.arrange(p0,p1, nrow = 1, top=textGrob("BNU1 Simulation",gp=gpar(fontsize=20,font=3)))

dev.off()


pdf(file="BNU1Stencil.pdf",
    width = 7,
    height = 7)
plot.GG.blocks.general.color.stencil(plot.data[["BNU1"]]$org.stencils[[paste(BNU1)]],
                                    "",
                                    Tree.full.list[["brain"]],
                                    stencil.plot.BNU1)
dev.off()

pdf(file="3mrmhsbm.pdf",
    width = 7,
    height = 7)
plot.GG.blocks.general.color.stencil(plot.data[["m3"]]$org.stencils[[paste(m3)]],
                                     "",
                                     Tree.full.list[["3m-7com"]],
                                     stencil.plot.m3)
dev.off()
