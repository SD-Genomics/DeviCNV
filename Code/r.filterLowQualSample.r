args<-commandArgs(TRUE)
batchTag<-args[1] 
inSampleInfoTxt<-args[2]
norReadDepthDir<-args[3] 
lqSampleDir<- args[4] 
MQList<-args[5]

library("ggplot2")

nameList <-read.table(paste0(inSampleInfoTxt), sep="\t", header=T)
sampleNameList<-nameList$Sample
sampleLen<-length(sampleNameList)

MQList<-unlist(strsplit(MQList,","))

for(MQ in MQList){
  dat<-read.table(paste0(norReadDepthDir,"/",batchTag,".readDepth.normalizedChrX.",MQ,".txt"), head=T, sep="\t", fill=T, quote = "", check.names=F)
  dat<-dat[dat$Chr!="Y",]

  out<-matrix(rep(NA,sampleLen*(sampleLen-1)),sampleLen-1,sampleLen)
  X<-c()
  Y<-c()
  sampleStartPos<-length(dat)-sampleLen+1
  sampleEndPos<-length(dat)
  for(i in seq(sampleStartPos,sampleEndPos)){
    subOut<-c()
    for(j in seq(sampleStartPos,sampleEndPos)){
      if(i!=j){
        x=as.numeric(as.character(unlist(dat[i])))
        y=as.numeric(as.character(unlist(dat[j])))
        corV=cor(x,y)
        X<-c(X,names(dat[i]))
        Y<-c(Y,corV)
        subOut<-c(subOut,corV)
      }
    }
    out[,i-sampleStartPos+1]<-subOut
  }
  
  colnames(out)<-sampleNameList
  outPercentile<-t(apply(out, 2, quantile, probes=c(0.25), na.rm=TRUE))
  lowQualLabel<-rep("HQ",sampleLen)
  lowQualLabel[outPercentile[,colnames(outPercentile)=="75%"]<0.7]<-"LQ"
  outPercentile<-cbind(rownames(outPercentile),outPercentile,lowQualLabel)
  colnames(outPercentile)<-c("Sample","0%","25%","50%","75%","100%","LowQualSample?")
  write.table(outPercentile, file=paste0(lqSampleDir,"/",batchTag,".All.lowQualitySampleTest.",MQ,".txt"), sep="\t", row.names=F, col.names=T,quote=F)
    
  datGG<-data.frame(x=X, y=Y)
  titleTXT<-paste0(batchTag,"\nCorrelation between samples (",MQ,")")
  ggplot(datGG,aes(X,Y))+geom_boxplot()+theme_bw()+scale_y_continuous(limits=c(-1,1))+labs(x="Samples",y="Correlation coefficient")+
    geom_hline(yintercept=0.7,linetype="dashed", size=0.3, color="red")+ggtitle(titleTXT)+theme(axis.text.x = element_text(angle = 90, hjust = 1))  
  ggsave(filename=paste0(lqSampleDir,"/",batchTag,".All.lowQualitySampleTest.",MQ,".pdf"), plot=last_plot(), width = (sampleLen*3+100), units= "mm")
}