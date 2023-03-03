rm(list=ls())
library(lme4)
library(randomForest)
library(caret)
library(epiR)
library(pROC)
library(geepack)
library(vioplot)
library(mgcv)
library(ggplot2)
library(itsadug)
library(tidytext)
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(ggpubr)
library(gridExtra)
library(lmerTest)
library(lqmm)
#function to display GA at sample
lpline=function(dat,tit){
 
  dat %>% 
    mutate(Group = Group2, 
           Main_Index= reorder_within(Main_Index,GA, Group)) %>% 
    ggplot() +
    geom_segment(aes(x=Main_Index, xend=Main_Index, y=GA, yend=Del_GA_Calc), color="gray")  + 
    geom_point(aes(x=Main_Index, y=GA,color="Amniocentesis",shape="Amniocentesis")) +
    geom_point(aes(x=Main_Index, y=Del_GA_Calc,color="Delivery",shape="Delivery")) +
    geom_hline(yintercept=c(24, 28, 32,37), col="blue") +
    labs(x="Patients",y="Gestational age (weeks)")+
    coord_flip()+
    theme_bw() +
    scale_x_reordered()+
    scale_y_continuous(breaks = seq(8,42,2),limits = c(6, 43))+
    scale_colour_manual(name="",values = c("grey38","darkgreen"),
                        labels=c("Gestational age at sample","Gestational age at delivery"))+
    scale_shape_manual(name="",values = c(16,6),
                       labels=c("Gestational age at sample","Gestational age at delivery"))+
    theme(legend.position="bottom")+
    facet_wrap(~Group, ncol=1, scale="free")+
    ggtitle(tit)
}

#load data

limss=c(37,37,34,34)
names(limss)<-c("PTL_34-37","PPROM_34-37","PTL_20-34","PPROM_20-34")

limss2=c(37,37,34,34)
names(limss2)<-c("PTL_NoInf","PPROM_NoInf","PTL_Inf","PPROM_Inf")


load("proteinandannotation2.RData")



pdf("Figure_GAsample.pdf")
lpline(ano,"")

pts=c(length(unique(ano$Main_Index[!duplicated(ano$Main_Index)&ano$Group2=="Control"])),
      length(unique(ano$Main_Index[!duplicated(ano$Main_Index)&ano$Group2=="PTL"])),
      length(unique(ano$Main_Index[!duplicated(ano$Main_Index)&ano$Group2=="PPROM"])))
smps=c(length((ano$Main_Index[ano$Group2=="Control"])),
       length((ano$Main_Index[ano$Group2=="PTL"])),length((ano$Main_Index[ano$Group2=="PPROM"])))
lbs=paste(c("Control","PTL","PPROM"),paste(c("N=","N=","N="),pts),
          paste(c("n=","n=","n="),smps),sep="\n")

vioplot(GA~Group2,ano,main="",names=lbs,xlab=NULL,ylab="Gestational age at sample (weeks)",
        col=c("blue","grey","red"))
dev.off()


ans=colnames(X)


#### changes with GA in controls lme models
#changes with GA overall
outga=NULL
pdf("cyt_TDvsGA_nolog.pdf")
#analyze diversity continous vs GA
for(an in ans){
  ano$Y=log2(X[,an])
  
  
  mod0=lmer(Y~1+(1|Main_Index),data=ano[ano$Group=="Control",],REML=FALSE)
  mod1=lmer(Y~GA+(1|Main_Index),data=ano[ano$Group=="Control",],REML=FALSE)
  #R2=performance::r2(mod1)$R2_marginal
  r=cor(ano[ano$Group=="Control","Y"],ano[ano$Group=="Control","GA"],method="spearman")
  outga=rbind(outga,c(summary(mod1)$coef["GA",c(1,5)],r))
  
  plme=anova(mod0,mod1)[["Pr(>Chisq)"]][2]
  plot(2^(Y)~GA,data=ano[ano$Group=="Control",],ylab=an,xlab="Gestational age (weeks)",
       main=paste("p=",round(plme,3)),log="y")
  mod1lme=mod1
  if(TRUE){
    mod1=bam(Y~s(GA, k=2)+ s(Main_Index, bs="re"), data=ano[ano$Group=="Control",], method="REML") 
    dat=ano[ano$Group=="Control",]
    ds=expand.grid(GA=seq(10,40,by=0.5)) 
    pred=get_predictions(mod1, cond=list(GA=ds$GA),rm.ranef=TRUE)
    pred$se=pred$CI/(1.96);pred$lwr=pred$fit-1.96*pred$se;pred$upr=pred$fit+1.96*pred$se;
    
    pred0=pred
    polygon(c(rev(pred0$GA), pred0$GA), 2^c(rev((pred0$lwr)), (pred0$upr)), col = 'grey', border = NA)
    points(2^(fit)~GA,data=pred,type="l",lwd=3,col="blue",
           ylab="",xlab="",xlim=c(10,40)) 
  }
  cat(an)
  cat("\n")
  
}
dev.off()
rownames(outga)<-ans
outga=data.frame(outga)
names(outga)[2]<-"P"
names(outga)[3]<-"r"
outga$q=p.adjust(outga$P,"fdr")
outga=outga[order(outga$P),]
write.csv(outga,file="cyt_vs_ga_r.csv")


### quantile regression models
#### controls GA quant reg
#changes with GA overall
outga=data.frame(GA=seq(10,40,1),Main_Index=0)
pdf("cyt_TDvsGA_qunatreg.pdf")
#analyze diversity continous vs GA
for(an in ans){
  ano$Y=log2(X[,an])
  
  fit.lqmm <- lqmm(fixed = Y~GA, random = ~ 1, group = Main_Index,	
                   data = ano[ano$Group=="Control",], tau = c(0.1,0.5,0.9))
  z=data.frame(predict(fit.lqmm,outga))
  names(z)=paste(an,c("10th","50th","95th"),sep="_")
  outga=cbind(outga,2^z)
  plot(outga$GA,2^z[,2],type="l",ylim=range(2^z),log="y",xlab="Gestational age (weeks)", ylab=paste("Concentration",an))
  points(outga$GA,2^z[,1],type="l",col="red")
  points(outga$GA,2^z[,3],type="l",col="red")
  
}
dev.off()
write.csv(outga,file="cyt_quantilesvsGA.csv")




# take only data collected preterm to compare among groups


ano=ano[ano$GA<=36.43,]
X=X[rownames(ano),]


pdf("Figure_cyt_longit.pdf")
#overall diveristy 
par(mfrow=c(2,2))
gps=levels(ano$Group2);cols=c("black","orange","purple");names(cols)<-gps

for(an in ans){
  ano$Y=log2(X[,an])
  #see if it was significant in Table S4
  gps=levels(ano$Group2);cols=c("black","orange","purple");names(cols)<-gps
  modnb=lmer(Y~Group2+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=ano,REML=FALSE)
  nbbm=summary(modnb)$coef[c("Group2PTL","Group2PPROM"),c("Pr(>|t|)")] 
  
  if(any(nbbm<c(0.052,0.0162))){#these are nominal p-values corresponding to q<0.1 in Table S4
    if(an!="VEGF"){
      plot(0,0,ylab=paste("log2",an,"concentration"),xlab="Gestational age (weeks)",
           main="",col=cols[as.character(ano$Group3)],ylim=c(quantile(ano$Y[ano$Group2=="PPROM"],0.3),quantile(ano$Y[ano$Group2=="PPROM"],0.9)),xlim=c(8,37))
    }else{
      plot(0,0,ylab=paste("log2",an,"concentration"),xlab="Gestational age (weeks)",
           main="",col=cols[as.character(ano$Group3)],ylim=c(quantile(ano$Y[ano$Group2=="PPROM"],0.1),quantile(ano$Y[ano$Group2=="PPROM"],0.9)),xlim=c(8,37))
      
    }

    for(g in gps){
      mod1=bam(Y~GA+ s(Main_Index, bs="re"), data=ano[ano$Group2==g,], method="REML") 
      ds=expand.grid(GA=seq(10,36.5,by=0.1),Main_Index="new");ds$GA2=ds$GA/40 
      pred=get_predictions(mod1, cond=list(GA=ds$GA),rm.ranef=TRUE)
      pred0=pred
      pred$se=pred$CI/(1.96);pred$lwr=pred$fit-1.96*pred$se;pred$upr=pred$fit+1.96*pred$se;
      points((fit)~GA,data=pred0[pred0$GA<ifelse(g%in%c("PTL_20-34","PPROM_20-34"),34,36.4),],type="l",lwd=3,col=cols[g],
             ylab="",xlab="",xlim=c(10,35))  
    }
    #interaction
    modnb=lmer(Y~Group2*GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=ano,REML=FALSE)
    nbb=summary(modnb)$coef[c("Group2PTL:GA2","Group2PPROM:GA2"),c("Pr(>|t|)")] 
    
    legend("topleft",lwd=rep(2,length(gps)),col=cols[gps],legend=paste(gps,c("",ifelse(nbb<0.05,"Interaction Sig.","Interaction NS"))),cex=0.9,ncol=1,bty="n")
    
  }
}

dev.off()



pdf("Figure_cyt_spagetti.pdf")

gps=levels(ano$Group2);cols=c("black","orange","purple");names(cols)<-gps

an="MCP1"
  
  ano$Y=log2(X[,an])
  #see if it was significant in Table S4
  gps=levels(ano$Group2);cols=c("darkgrey","orange","purple");names(cols)<-gps
  modnb=lmer(Y~Group2+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=ano,REML=FALSE)
  nbbm=summary(modnb)$coef[c("Group2PTL","Group2PPROM"),c("Pr(>|t|)")] 
  
  if(any(nbbm<c(0.052,0.0162))){
    plot(0,0,ylab=paste("log2",an,"concentration"),xlab="Gestational age (weeks)",
         main="",ylim=c(quantile(ano$Y,0.01),quantile(ano$Y,0.99)),xlim=c(8,37))
  
    for(p in unique(ano$Main_Index)){
      tm=ano[ano$Main_Index==p,]    
      points(tm$GA,tm$Y,type="l", col=cols[as.character(tm$Group2)])  
    }
    
    
  }


dev.off()



pdf("cytokines_by_groups.pdf")
#overall diveristy 
for(an in ans){
  ano$Y=log2(X[,an])
  ano$Y2=X[,an]
  ano[,an]<-X[,an]
  
  stat.test <- compare_means(
    Y2 ~ Group, data = ano,
    method = "t.test"
  )%>%mutate(y.position=log2(max(ano[,an])*1.5))
  
  mod1=lmer(Y~Group+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=ano,REML=FALSE)
  mod0=lmer(Y~GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=ano,REML=FALSE)
  ps=round(anova(mod0,mod1)[["Pr(>Chisq)"]][2],4)
  
  stat.test$p=ps
  
  p1 <- ggviolin(
    ano, x = "Group", y = an,
    color = "Group", add = "jitter", xlab = paste(""), ylab = an)+
    stat_pvalue_manual(data = stat.test,label = "p")+ scale_color_manual(values=c("blue", "orange"))+
    theme(legend.position="none",axis.text.x=element_text(face = "bold", size = 8))+  
    scale_y_continuous(trans='log2')+ 
    stat_summary(fun = "median",
                 geom = "crossbar", 
                 width = 0.5,
                 colour = c("black"),linetype="longdash")
  
  
  #analyze diversity binary
  stat.test <- compare_means(
    Y ~ Group2, data = ano,
    method = "t.test"
  )
  stat.test=stat.test[-3,]
  stat.test$y.position=max(ano$Y)*1.5+(0:1)
  ps=NULL
  for(gps in levels(ano$Group2)[-1]){
    sdat=ano[ano$Group2%in%c("Control",gps),]
    mod1=lmer(Y~Group2+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=sdat,REML=FALSE)
    mod0=lmer(Y~GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=sdat,REML=FALSE)
    ps=c(ps,round(anova(mod0,mod1)[["Pr(>Chisq)"]][2],4))
  }
  stat.test$p=ps
  
  p2 <- ggviolin(
    ano, x = "Group2", y = an,
    color = "Group2", add = "jitter", xlab = paste(""), ylab = an)+
    stat_pvalue_manual(data = stat.test,label = "p")+ scale_color_manual(values=c("blue","darkgrey","red"))+
    theme(legend.position="none",axis.text.x=element_text(face = "bold", size = 8))+  
    scale_y_continuous(trans='log2')+ 
    stat_summary(fun = "median",
                 geom = "crossbar", 
                 width = 0.5,
                 colour = c("black"),linetype="longdash")
  
  #by PTL PPROM and IR
  
  stat.test <- compare_means(
    Y ~ Group3, data = ano,
    method = "t.test"
  )
  stat.test=stat.test[-c(5:10),]
  stat.test$y.position=max(ano$Y)*1.5+(0:3)
  
  ps=NULL;out3g=NULL
  for(gps in levels(ano$Group3)[-1]){
    sdat=ano[ano$Group3%in%c("Control",gps),]
    sdat=sdat[sdat$GA<=limss[gps],]
    
    mod1=lmer(Y~Group3+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=sdat,REML=FALSE)
    mod0=lmer(Y~GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=sdat,REML=FALSE)
    ps=c(ps,round(anova(mod0,mod1)[["Pr(>Chisq)"]][2],4))
  }
  stat.test$p=ps
  
  
  p3 <- ggviolin(
    ano, x = "Group3", y = an,
    color = "Group3", add = "jitter", xlab = paste(""), ylab = an)+
    stat_pvalue_manual(data = stat.test,label = "p")+ scale_color_manual(values=c("blue","grey44","grey","red4","red"))+
    theme(legend.position="none",axis.text.x=element_text(face = "bold", size = 8))+  
    scale_y_continuous(trans='log2')+ 
    stat_summary(fun = "median",
                 geom = "crossbar", 
                 width = 0.5,
                 colour = c("black"),linetype="longdash")
  layout <- rbind(c(1, 2),
                  c(3, 3))
  grid.arrange(p1,p2,p3, layout_matrix=layout)
  
}

dev.off()




##
#overall differences between sPTB and COntrols
NN=dim(ano)[1]
set.seed(1)
res=NULL
for(gps in levels(ano$Group)[-1]){
  sdat=ano[ano$Group%in%c("Control",gps),]
  
  tres=NULL 
  for (an in colnames(X)){ #for each protein
    #NB binary
    sdat$Y=log2(X[rownames(sdat),an])
    
    modnb=lmer(Y~Group+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=sdat,REML=FALSE)
    nbb=summary(modnb)$coef[2,c("Estimate","Pr(>|t|)")]
    tres=rbind(tres,nbb)
  }
  tres=cbind(tres,p.adjust(tres[,2],"fdr"))
  colnames(tres)<-c(paste("Coef_",gps,"vsControl",sep=""),paste("P_",gps,"vsControl",sep=""),paste("Q_",gps,"vsControl",sep=""))
  res=cbind(res,tres)
}
res1=data.frame(ID=colnames(X),res)

res=NULL
for(gps in levels(ano$Group)[-1]){
  sdat=ano[ano$Group%in%c("Control",gps),]
  
  tres=NULL 
  for (an in colnames(X)){ #for each microbial species
    #NB binary
    sdat$Y=log2(X[rownames(sdat),an])
    
    modnb=lmer(Y~Group+GA2+(1|Main_Index),data=sdat,REML=FALSE)
    nbb=summary(modnb)$coef[2,c("Estimate","Pr(>|t|)")]
    tres=rbind(tres,nbb)
  }
  tres=cbind(tres,p.adjust(tres[,2],"fdr"))
  colnames(tres)<-c(paste("Coef_",gps,"vsControl",sep=""),paste("P_",gps,"vsControl",sep=""),paste("Q_",gps,"vsControl",sep=""))
  res=cbind(res,tres)
}
res2=data.frame(ID=colnames(X),res)

res=cbind(res1,res2)
res=res[order(res[,3]),]
write.csv(res,file="cyt_results_overall_withandwithoutadj.csv")



#differences subtypes of sPTB (PPROM; PTL)
NN=dim(ano)[1]
set.seed(1)
res=NULL
for(gps in levels(ano$Group2)[-1]){
  sdat=ano[ano$Group2%in%c("Control",gps),]
  tres=NULL 
  for (an in colnames(X)){ #for each microbial species
    #NB binary
    sdat$Y=log2(X[rownames(sdat),an])
    modnb=lmer(Y~Group2+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=sdat,REML=FALSE)
    nbb=summary(modnb)$coef[2,c("Estimate","Pr(>|t|)")] 
    tres=rbind(tres,nbb)
  }
  tres=cbind(tres,p.adjust(tres[,2],"fdr"))
  colnames(tres)<-c(paste("Coef_",gps,"vsControl",sep=""),paste("P_",gps,"vsControl",sep=""),paste("Q_",gps,"vsControl",sep=""))
  res=cbind(res,tres)
}
res=data.frame(ID=colnames(X),res)
res=res[order(res$P_PTL),]

write.csv(res,file="cyrres_by_phenotype.csv")



#Changes by  early and late PPROM and PTL
NN=dim(ano)[1]
set.seed(1)
res=NULL
for(gps in levels(ano$Group3)[-1]){
  sdat=ano[ano$Group3%in%c("Control",gps),]
  sdat=sdat[sdat$GA<=limss[gps],]
  tres=NULL 
  for (an in colnames(X)){ #for each protein
    #NB binary
    sdat$Y=log2(X[rownames(sdat),an])
    modnb=lmer(Y~Group3+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index),data=sdat,REML=FALSE)
    nbb=summary(modnb)$coef[2,c("Estimate","Pr(>|t|)")] 
    tres=rbind(tres,nbb)
  }
  tres=cbind(tres,p.adjust(tres[,2],"fdr"))
  colnames(tres)<-c(paste("Coef_",gps,"vsControl",sep=""),paste("P_",gps,"vsControl",sep=""),paste("Q_",gps,"vsControl",sep=""))
  res=cbind(res,tres)
}
res=data.frame(ID=colnames(X),res)

write.csv(res,file="cyrres_by_phenotypeGAD.csv")



#GA del with PPROM
pdf("CytvsGADel.pdf")
#"max","last","slope"
par(mfrow=c(2,2))
for(outcom in c("PTL","PPROM")){
  for(cutm in c(28)){
    
    for(suma in c("last")){
      ano2=ano[ano$GA<=cutm&ano$Del_GA_Calc>cutm,]
      if(outcom=="PTL"){ano2=ano2[ano2$Group2%in%c("PTL"),]}
      if(outcom=="PPROM"){ano2=ano2[ano2$Group2%in%c("PPROM"),]}
      #ano2=ano2[!(ano2$Group=="Control"&ano2$Lab_Sp==0),]
      #
      if(suma=="slope"){
        ano2=ano2[ano2$Main_Index%in%names(table(ano2$Main_Index))[table(ano2$Main_Index)>=2],]
      }
      
      if(suma=="last"){
        ano2=ano2[order(ano2$GA,decreasing =TRUE),]
        ano2=ano2[!duplicated(ano2$Main_Index),]
      }
      
      ano2$Main_Index=factor(as.character(ano2$Main_Index))
      eset=t(X[rownames(ano2),])
      
      if(suma%in%c("mean","last")){
        eset=t(apply(eset,1,function(x){tapply(x,ano2$Main_Index,mean)}))
      }
      if(suma=="max"){
        eset=t(apply(eset,1,function(x){tapply(x,ano2$Main_Index,max)}))
      }
      if(suma=="slope"){
        esets=matrix(NA,dim(eset)[1],length(levels(ano2$Main_Index)))
        rownames(esets)<-rownames(eset)
        colnames(esets)<-levels(ano2$Main_Index)
        for(pat in levels(ano2$Main_Index)){
          tmp=ano2[ano2$Main_Index==pat,]
          for(as in rownames(eset)){
            tmp$Y <-eset[as,rownames(tmp)]
            esets[as,pat]<-lm(Y~GA,tmp)$coef[2]
          }
        }
        eset=esets
      }
      
      
      ano2=ano2[!duplicated(ano2$Main_Index),]
      ano2=ano2[match(colnames(eset),ano2$Main_Index),]
      
      for (div in ans){
        mydiv=eset[div,]
        GAD=ano2$Del_GA_Calc
        mod=lm(mydiv~GAD)
        plot(GAD,mydiv,xlab="Gestational age at delivery",ylab=div,pch=19,cex.main=0.8,main=paste("Ouctome: ", outcom,"; Samples 8-", cutm," wks; Summary: ",suma,sep=""))
        abline(mod$coef,lwd=2)
        r=round(cor(mydiv,GAD,method="spearman"),2)
        legend("topleft",c(paste("R=",r,sep=""),paste("p=",round(cor.test(GAD,mydiv,method="spearman")$p.value,3),sep="")))
        
      }
    }
  }
}
dev.off()




####predictions
pdf(paste("cyt_ROCs_last_34_matchar.pdf",sep="_"))
#"max","last","slope"
par(mfrow=c(2,2))
for(outcom in c("PTL_PPROM","PTL","PPROM","PPROM<34","PTL<34")){
  for(cutm in c(28)){
    for(suma in c("last")){
      ano2=ano[ano$GA<=cutm&ano$Del_GA_Calc>cutm,]
      if(outcom=="PTL"){ano2=ano2[ano2$Group2%in%c("Control","PTL"),]}
      if(outcom=="PPROM"){ano2=ano2[ano2$Group2%in%c("Control","PPROM"),]}
      if(outcom=="PPROM<34"){ano2=ano2[ano2$Group2=="Control"|(ano2$Group2=="PPROM"&ano2$Del_GA_Calc<34),]}
      if(outcom=="PTL<34"){ano2=ano2[ano2$Group2=="Control"|(ano2$Group2=="PTL"&ano2$Del_GA_Calc<34),]}
      
      #
      if(suma=="slope"){
        ano2=ano2[ano2$Main_Index%in%names(table(ano2$Main_Index))[table(ano2$Main_Index)>=2],]
      }
      
      if(suma=="last"){
        ano2=ano2[order(ano2$GA,decreasing =TRUE),]
        ano2=ano2[!duplicated(ano2$Main_Index),]
      }
      
      ano2$Main_Index=factor(as.character(ano2$Main_Index))
      eset=t(X[rownames(ano2),])
      
      if(suma%in%c("mean","last")){
        eset=t(apply(eset,1,function(x){tapply(x,ano2$Main_Index,mean)}))
      }
      if(suma=="max"){
        eset=t(apply(eset,1,function(x){tapply(x,ano2$Main_Index,max)}))
      }
      if(suma=="slope"){
        esets=matrix(NA,dim(eset)[1],length(levels(ano2$Main_Index)))
        rownames(esets)<-rownames(eset)
        colnames(esets)<-levels(ano2$Main_Index)
        for(pat in levels(ano2$Main_Index)){
          tmp=ano2[ano2$Main_Index==pat,]
          for(as in rownames(eset)){
            tmp$Y <-eset[as,rownames(tmp)]
            esets[as,pat]<-lm(Y~GA,tmp)$coef[2]
          }
        }
        eset=esets
        
      }
      
      
      ano2=ano2[!duplicated(ano2$Main_Index),]
      ano2=ano2[match(colnames(eset),ano2$Main_Index),]
      ano2$G=ifelse(ano2$Group=="PTL_PPROM","D","C")
      
      
      #models with protein data
      set.seed(200)
      perfs=NULL
      cx=NULL
      outcs<-predps<-NULL
      parts=createMultiFolds(ano2$G, k=10,times = 1)
      for(j in 1:10){
        trainIndex <-parts[[j]] 
        train=data.frame(G=factor(ano2$G[trainIndex]),t(eset[,trainIndex]))
        test=data.frame(G=factor(ano2$G[-trainIndex]),t(eset[,-trainIndex]))
        #m2=ML_RF(train,test,nfs=c(20,50),cut=0.9)
        modb= randomForest(train[,-1], y=train$G,   ntree=1000,importance=FALSE)
        out=predict(modb,test,type="prob")
        cut=quantile(out[test$G=="C","D"],0.9)
        pred=predict(modb,test,type="prob")[,"D"]
        Yh=ifelse(pred>cut,1,0)
        a=table(Yh,test$G)
        r=epi.tests(a[c(2,1),c(2,1)])
        outcs=c(outcs,as.numeric(test$G=="D"))
        predps=c(predps,pred)
        perfs=c(perfs,auc(roc(response=as.numeric(test$G=="D"),predictor=pred)))
      }
      
      aucs=round(ci.auc(roc(response=outcs,predictor=predps)),3)
      
      RC=roc(response=outcs,predictor=predps)
      plot(1-RC$specificities,RC$sensitivities,col="black",type="l",lwd=2,main=outcom,xlab="False positive rate",
           ylab="Sensitivity")
      abline(0,1)
      
      #maternal characteristics alone
      esetm=rbind(t(as.matrix(ano2[,c("BMI60","Parity","Preterm2","Age")])))
      
      set.seed(200)
      perfs=NULL
      cx=NULL
      outcs<-predps<-NULL
      parts=createMultiFolds(ano2$G, k=10,times = 1)
      for(j in 1:10){
        trainIndex <-parts[[j]] 
        train=data.frame(G=factor(ano2$G[trainIndex]),t(esetm[,trainIndex]))
        test=data.frame(G=factor(ano2$G[-trainIndex]),t(esetm[,-trainIndex]))
        #m2=ML_RF(train,test,nfs=c(20,50),cut=0.9)
        modb= randomForest(train[,-1], y=train$G,   ntree=1000,importance=FALSE)
        out=predict(modb,test,type="prob")
        cut=quantile(out[test$G=="C","D"],0.9)
        pred=predict(modb,test,type="prob")[,"D"]
        Yh=ifelse(pred>cut,1,0)
        a=table(Yh,test$G)
        r=epi.tests(a[c(2,1),c(2,1)])
        outcs=c(outcs,as.numeric(test$G=="D"))
        predps=c(predps,pred)
        perfs=c(perfs,auc(roc(response=as.numeric(test$G=="D"),predictor=pred)))
      }
      
      
      aucs2=round(ci.auc(roc(response=outcs,predictor=predps)),3)
      RC1=roc(response=outcs,predictor=predps)
      points(1-RC1$specificities,RC1$sensitivities,col="blue",type="l",lwd=2,main=outcom,xlab="False positive rate",
             ylab="Sensitivity")
      
      #proteins and maternal characteristics
      esetm=rbind(eset,t(as.matrix(ano2[,c("BMI60","Parity","Preterm2","Age")])))
      
      
      set.seed(200)
      perfs=NULL
      outcs<-predps<-NULL
      parts=createMultiFolds(ano2$G, k=10,times = 1)
      for(j in 1:10){
        trainIndex <-parts[[j]] 
        train=data.frame(G=factor(ano2$G[trainIndex]),t(esetm[,trainIndex]))
        test=data.frame(G=factor(ano2$G[-trainIndex]),t(esetm[,-trainIndex]))
        
        modb= randomForest(train[,-1], y=train$G,   ntree=1000,importance=FALSE)
        out=predict(modb,test,type="prob")
        cut=quantile(out[test$G=="C","D"],0.9)
        pred=predict(modb,test,type="prob")[,"D"]
        Yh=ifelse(pred>cut,1,0)
        a=table(Yh,test$G)
        r=epi.tests(a[c(2,1),c(2,1)])
        outcs=c(outcs,as.numeric(test$G=="D"))
        predps=c(predps,pred)
        #ci.auc(roc(response=as.numeric(test$G=="D"),predictor=pred))
        perfs=c(perfs,auc(roc(response=as.numeric(test$G=="D"),predictor=pred)))
      }
      
      aucs3=round(ci.auc(roc(response=outcs,predictor=predps)),3)
      RC3=roc(response=outcs,predictor=predps)
      points(1-RC3$specificities,RC3$sensitivities,col="red",type="l",lwd=2,main=outcom,xlab="False positive rate",
             ylab="Sensitivity")
      
      legs=c(paste("AUC=",aucs[2],"(",aucs[1],"-",aucs[3],")",sep=""),paste("AUC=",aucs2[2],"(",aucs2[1],"-",aucs2[3],")",sep=""),paste("AUC=",aucs3[2],"(",aucs3[1],"-",aucs3[3],")",sep=""))
      legend("bottomright",legs,lwd=c(2,2),col=c("black","blue","red"),cex=0.7,bty ="n")
      legend("topleft",bty ="n",legend=paste("p=",round(roc.test(RC3,RC1,method = "delong")$p.v,3)))
    }
  }
}
dev.off()




##30 weeks cut-off


####predictions
pdf(paste("cyt_ROCs_last_30_matchar.pdf",sep="_"))
#"max","last","slope"
par(mfrow=c(2,2))
for(outcom in c("PTL_PPROM","PTL","PPROM","PPROM<30","PTL<30")){
  for(cutm in c(24)){
    
    for(suma in c("last")){
      ano2=ano[ano$GA<=cutm&ano$Del_GA_Calc>cutm,]
      if(outcom=="PTL"){ano2=ano2[ano2$Group2%in%c("Control","PTL"),]}
      if(outcom=="PPROM"){ano2=ano2[ano2$Group2%in%c("Control","PPROM"),]}
      if(outcom=="PPROM<30"){ano2=ano2[ano2$Group2=="Control"|(ano2$Group2=="PPROM"&ano2$Del_GA_Calc<30),]}
      if(outcom=="PTL<30"){ano2=ano2[ano2$Group2=="Control"|(ano2$Group2=="PTL"&ano2$Del_GA_Calc<30),]}
      
      #
      if(suma=="slope"){
        ano2=ano2[ano2$Main_Index%in%names(table(ano2$Main_Index))[table(ano2$Main_Index)>=2],]
      }
      
      if(suma=="last"){
        ano2=ano2[order(ano2$GA,decreasing =TRUE),]
        ano2=ano2[!duplicated(ano2$Main_Index),]
      }
      
      ano2$Main_Index=factor(as.character(ano2$Main_Index))
      eset=t(X[rownames(ano2),])
      
      if(suma%in%c("mean","last")){
        eset=t(apply(eset,1,function(x){tapply(x,ano2$Main_Index,mean)}))
      }
      if(suma=="max"){
        eset=t(apply(eset,1,function(x){tapply(x,ano2$Main_Index,max)}))
      }
      if(suma=="slope"){
        esets=matrix(NA,dim(eset)[1],length(levels(ano2$Main_Index)))
        rownames(esets)<-rownames(eset)
        colnames(esets)<-levels(ano2$Main_Index)
        for(pat in levels(ano2$Main_Index)){
          tmp=ano2[ano2$Main_Index==pat,]
          for(as in rownames(eset)){
            tmp$Y <-eset[as,rownames(tmp)]
            esets[as,pat]<-lm(Y~GA,tmp)$coef[2]
          }
        }
        eset=esets
        
      }
      
      ano2=ano2[!duplicated(ano2$Main_Index),]
      ano2=ano2[match(colnames(eset),ano2$Main_Index),]
      ano2$G=ifelse(ano2$Group=="PTL_PPROM","D","C")
      
      
      #proteins alone
      set.seed(200)
      perfs=NULL
      cx=NULL
      outcs<-predps<-NULL
      parts=createMultiFolds(ano2$G, k=10,times = 1)
      for(j in 1:10){
        trainIndex <-parts[[j]] 
        train=data.frame(G=factor(ano2$G[trainIndex]),t(eset[,trainIndex]))
        test=data.frame(G=factor(ano2$G[-trainIndex]),t(eset[,-trainIndex]))
        #m2=ML_RF(train,test,nfs=c(20,50),cut=0.9)
        modb= randomForest(train[,-1], y=train$G,   ntree=1000,importance=FALSE)
        out=predict(modb,test,type="prob")
        cut=quantile(out[test$G=="C","D"],0.9)
        pred=predict(modb,test,type="prob")[,"D"]
        Yh=ifelse(pred>cut,1,0)
        a=table(Yh,test$G)
        r=epi.tests(a[c(2,1),c(2,1)])
        outcs=c(outcs,as.numeric(test$G=="D"))
        predps=c(predps,pred)
        #ci.auc(roc(response=as.numeric(test$G=="D"),predictor=pred))
        perfs=c(perfs,auc(roc(response=as.numeric(test$G=="D"),predictor=pred)))
      }
      
      aucs=round(ci.auc(roc(response=outcs,predictor=predps)),3)
      
      RC=roc(response=outcs,predictor=predps)
      plot(1-RC$specificities,RC$sensitivities,col="black",type="l",lwd=2,main=outcom,xlab="False positive rate",
           ylab="Sensitivity")
      abline(0,1)
      
      #maternal characteristics alone
      esetm=rbind(t(as.matrix(ano2[,c("BMI60","Parity","Preterm2","Age")])))
      
      set.seed(200)
      perfs=NULL
      cx=NULL
      outcs<-predps<-NULL
      parts=createMultiFolds(ano2$G, k=10,times = 1)
      for(j in 1:10){
        trainIndex <-parts[[j]] 
        train=data.frame(G=factor(ano2$G[trainIndex]),t(esetm[,trainIndex]))
        test=data.frame(G=factor(ano2$G[-trainIndex]),t(esetm[,-trainIndex]))
        #m2=ML_RF(train,test,nfs=c(20,50),cut=0.9)
        modb= randomForest(train[,-1], y=train$G,   ntree=1000,importance=FALSE)
        out=predict(modb,test,type="prob")
        cut=quantile(out[test$G=="C","D"],0.9)
        pred=predict(modb,test,type="prob")[,"D"]
        Yh=ifelse(pred>cut,1,0)
        a=table(Yh,test$G)
        r=epi.tests(a[c(2,1),c(2,1)])
        outcs=c(outcs,as.numeric(test$G=="D"))
        predps=c(predps,pred)
        #ci.auc(roc(response=as.numeric(test$G=="D"),predictor=pred))
        perfs=c(perfs,auc(roc(response=as.numeric(test$G=="D"),predictor=pred)))
      }
      
      
      aucs2=round(ci.auc(roc(response=outcs,predictor=predps)),3)
      RC1=roc(response=outcs,predictor=predps)
      points(1-RC1$specificities,RC1$sensitivities,col="blue",type="l",lwd=2,main=outcom,xlab="False positive rate",
             ylab="Sensitivity")
      
      #proteins with matchar
      esetm=rbind(eset,t(as.matrix(ano2[,c("BMI60","Parity","Preterm2","Age")])))
      
      
      set.seed(200)
      perfs=NULL
      outcs<-predps<-NULL
      parts=createMultiFolds(ano2$G, k=10,times = 1)
      for(j in 1:10){
        trainIndex <-parts[[j]] 
        train=data.frame(G=factor(ano2$G[trainIndex]),t(esetm[,trainIndex]))
        test=data.frame(G=factor(ano2$G[-trainIndex]),t(esetm[,-trainIndex]))
        
        #m2=ML_RF(train,test,nfs=c(20,50),cut=0.9)
        modb= randomForest(train[,-1], y=train$G,   ntree=1000,importance=FALSE)
        out=predict(modb,test,type="prob")
        cut=quantile(out[test$G=="C","D"],0.9)
        pred=predict(modb,test,type="prob")[,"D"]
        Yh=ifelse(pred>cut,1,0)
        a=table(Yh,test$G)
        r=epi.tests(a[c(2,1),c(2,1)])
        outcs=c(outcs,as.numeric(test$G=="D"))
        predps=c(predps,pred)
        #ci.auc(roc(response=as.numeric(test$G=="D"),predictor=pred))
        perfs=c(perfs,auc(roc(response=as.numeric(test$G=="D"),predictor=pred)))
      }
      
      aucs3=round(ci.auc(roc(response=outcs,predictor=predps)),3)
      RC3=roc(response=outcs,predictor=predps)
      points(1-RC3$specificities,RC3$sensitivities,col="red",type="l",lwd=2,main=outcom,xlab="False positive rate",
             ylab="Sensitivity")
      legs=c(paste("AUC=",aucs[2],"(",aucs[1],"-",aucs[3],")",sep=""),paste("AUC=",aucs2[2],"(",aucs2[1],"-",aucs2[3],")",sep=""),paste("AUC=",aucs3[2],"(",aucs3[1],"-",aucs3[3],")",sep=""))
      legend("bottomright",legs,lwd=c(2,2),col=c("black","blue","red"),cex=0.7,bty ="n")
      legend("topleft",bty ="n",legend=paste("p=",round(roc.test(RC3,RC1,method = "delong")$p.v,3)))
    }
  }
}
dev.off()




