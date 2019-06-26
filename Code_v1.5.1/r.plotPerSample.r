args<-commandArgs(TRUE)
batchTag<-args[1] 
inSample<-args[2]
RDRatioDir<-args[3] 
plotdir<- args[4] 
CBSdir<-args[5]
PoolList<-args[6]
MQList<-args[7]
dupdelList<-args[8]

RvalTh<-0.8
rangeX<-1 ##How many gene you want to plot in one figure.

library("ggplot2")
library("PSCBS")

PoolList<-unlist(strsplit(PoolList,","))
MQList<-unlist(strsplit(MQList,","))
dupdelList<-unlist(strsplit(dupdelList,","))

for(MQ in MQList){
  for(dupdelTh in dupdelList){
    dupTh<-as.numeric(unlist(strsplit(dupdelTh,"_"))[1])
    delTh<-as.numeric(unlist(strsplit(dupdelTh,"_"))[2])
    
    inFile=paste0(RDRatioDir,batchTag,".readDepthRatioFromLRModel.",MQ,".dupdelTh",dupdelTh,".txt")
    RDTable<-read.table(inFile, head=T, sep="\t", fill=T, quote = "", check.names=F)
    RDTable<-RDTable[order(RDTable$Amplicon_Start),]
    sortedRDTable<-c()
    for(chrom in c(1:22, "X", "Y", "MT")){
      sortedRDTable<-rbind(sortedRDTable,RDTable[RDTable$Chr==chrom,])
    }
    RDTable<-sortedRDTable
    
    AllPlot=TRUE
    if(AllPlot==TRUE){
      ### All #####################
      geneList<-unique(RDTable$Gene)
      geneLen<-NROW(geneList)
      
      medCov<-unique(RDTable[RDTable$Type=="MedianRD",][c("Pool",inSample)])
      colnames(medCov)<-c("Pool","MedianRD")
      
      ampLen<-length(unique(RDTable$Amplicon_ID))
      xvalues<-rep(NA, ampLen)
      CNMvalues<-rep(NA, ampLen)
      pvalues<-rep(NA, ampLen)
      rvalues<-rep(NA, ampLen)
      cnvtypes<-rep(NA, ampLen)
      poolvalues<-rep(19, ampLen)
      alphavalues<-rep(1.0, ampLen)
      sizevalues<-rep(0.8, ampLen)
      genomicPos<-rep(NA, ampLen)
      chrvalues<-rep(NA, ampLen)
      chrlabels<-rep(NA, ampLen)
      gaps<-c()
      
      xstart<-1
      start<-1
      for(x in 1:geneLen){
        gene=geneList[x]
        RDTable.gene<-RDTable[RDTable$Gene==gene,]
        ampLen.gene<-length(unique(RDTable.gene$Amplicon_ID))
        end<-start+ampLen.gene-1
        xend<-xstart+10*(ampLen.gene-1)
        
        xvalues[start:end]<-seq(from=xstart, to=xend, by=10)
        CNMvalues[start:end]<-as.numeric(as.matrix(RDTable.gene[RDTable.gene$Type=="CN_M",inSample]))
        pvalues[start:end]<-as.numeric(as.matrix(RDTable.gene[RDTable.gene$Type=="Pvalue",inSample]))
        rvalues[start:end]<-as.numeric(as.matrix(RDTable.gene[RDTable.gene$Type=="RegRvalue",inSample]))
        cnvtypes[start:end]<-as.character(as.matrix(RDTable.gene[RDTable.gene$Type=="CNVType",inSample]))
        poolvalues[start:end]<-as.character(as.matrix(RDTable.gene[RDTable.gene$Type=="CN_M","Pool"]))  
        chrvalues[start:end]<-as.character(RDTable.gene[RDTable.gene$Type=="CN_M","Chr"])
        chrlabels[start:end]<-as.character(RDTable.gene[RDTable.gene$Type=="CN_M","Chr"])
        
        genomicPos[start:end]<-as.numeric(RDTable.gene[RDTable.gene$Type=="CN_M","Amplicon_Start"])
        gaps.chr<-chrvalues[start]
        gaps.chr<-gsub("X",23,gaps.chr)
        gaps.chr<-gsub("Y",24,gaps.chr)
        gaps.chr<-gsub("MT",25,gaps.chr)
        gaps.chr<-as.numeric(gaps.chr)
        gaps.start<-min(genomicPos[start:end])
        gaps.end<-max(genomicPos[start:end])+1
        gaps.length<-(gaps.end-gaps.start+1)
        gaps<-rbind(gaps, c(gaps.chr,gaps.start,gaps.end,gaps.length))
        
        xstart<-xend+1
        start<-end+1
        x<-x+1
      }
      
      if(length(pvalues[is.na(pvalues)])==length(pvalues)){
        print("All p-values are nan!!! impossible to draw plot")
      }else{
        
        MAXV=3
        MINV=(-2)

	chrvalues[chrvalues=="X"]<-23
        chrvalues[chrvalues=="Y"]<-24
        chrvalues[chrvalues=="MT"]<-25

        copynumber<-CNMvalues
	copynumber[copynumber>2^(MAXV)]<-2^(MAXV)*1.0
        copynumber[cnvtypes=="faultyAmp"]<-NA
        copynumber[cnvtypes=="faultySample"]<-NA
        
        CBSDat<-data.frame(chromosome=as.numeric(chrvalues), x=genomicPos, y=copynumber)
        CBSDat<-CBSDat[!is.na(CBSDat$y),]
        
        colnames(gaps)<-c("chromosome","start","end","length")	
        gaps2<-data.frame(gaps)
        if(NROW(gaps2)>1){
          for(i in seq(NROW(gaps2),2)){
            end1<-gaps2$end[i-1]
            start2<-gaps2$start[i]
            if(gaps2$chromosome[i-1]==gaps2$chromosome[i]){
              if(end1>start2){
                gaps2[i-1,3]<-gaps2$end[i]
                gaps2[i-1,4]<-gaps2[i-1,3]-gaps2[i-1,2]+1
                gaps2<-gaps2[-i,]
              }
            }
          }
          knownSegments<-gapsToSegments(gaps2)
          fit <- segmentByCBS(CBSDat,knownSegments=knownSegments)
        }else{
          fit <- segmentByCBS(CBSDat)
        }
        CBSFileName=paste0(CBSdir,"/",inSample,".CBS.",MQ,".dupdelTh",dupdelTh)
        CBSFile=paste0(CBSdir,"/",inSample,".CBS.",MQ,".dupdelTh",dupdelTh,".tsv")
        if(file.exists(CBSFile)){
          file.remove(CBSFile)
        }
        pathname<-writeSegments(fit, name=CBSFileName, simplify=TRUE)


        ###outlier handlin######
        temp.poolvalues<-poolvalues
        
        temp.CN_M<-CNMvalues
        temp.CN_M<-log2(temp.CN_M+0.000000001)
        temp.CN_M[temp.CN_M>MAXV]<-MAXV*1.0
        temp.CN_M[temp.CN_M<MINV]<-MINV*1.0
        CNMvalues<-temp.CN_M    
         
        poolInfo<-c()
        refCexList<-c(15,16,17,18)
        poolCexList<-c()
        i<-1
        for(Pool in PoolList){
          poolCheck<-(medCov$Pool %in% c(Pool))
          if(length(poolCheck[poolCheck==TRUE])>0){
            if(as.numeric(as.matrix(medCov[medCov$Pool==Pool,2]))<50){
              temp.poolvalues[temp.poolvalues==Pool]<-"LowMedRD"
            }
          }
          poolInfo<-cbind(poolInfo, paste0(Pool,":",medCov[medCov$Pool==Pool,2]))
          poolCexList<-cbind(poolCexList,refCexList[i])
          i<-i+1
        }
        
        chr.now=""
        for(chrPos in c(1:length(chrlabels))){
          if(is.na(chrlabels[chrPos])){
            chrlabels[chrPos]<-""
          }
          if(chr.now==chrlabels[chrPos]){
            chrlabels[chrPos]<-""
          }
          if(chrlabels[chrPos]!=""){
            chr.now<-chrlabels[chrPos]
          }
        }

        temp.poolvalues[rvalues<RvalTh]<-paste0("LowRval(<",RvalTh,")")
        temp.poolvalues[is.na(pvalues)]<-"Faulty"
        poolvalues<-factor(temp.poolvalues,levels=c(PoolList,"Faulty","LowMedRD",paste0("LowRval(<",RvalTh,")")))

	alphavalues[poolvalues=="Faulty"]<-0.2
        alphavalues[poolvalues=="LowMedRD"]<-0.6
        alphavalues[poolvalues=="LowRval"]<-0.6
        dat<-data.frame(xvalues, CNMvalues, pvalues, poolvalues, alphavalues, sizevalues, chrlabels)
        
        medCovTXT=paste0(poolInfo,collapse=", ")
        titleTXT=paste0(inSample," (",MQ,")\n",medCovTXT)
        ggplot(dat, aes(x=dat$xvalues,y=dat$CNMvalues, shape=dat$poolvalues, color=dat$pvalues, alpha=dat$alphavalues))+
          scale_x_continuous(breaks=c())+ 
          scale_y_continuous(limits=c(-3,3), breaks=seq(MINV,MAXV,1))+ 
          labs(x="All Genes", y="Log2(observed read depth/expected read depth)")+theme_bw()+	 
          ggtitle(titleTXT)+
          theme(plot.title = element_text(size=20))+
          geom_point(size=0.8)+scale_alpha(guide = 'none', limits=c(0,1))+
          scale_shape_manual("Pool",limits=c(PoolList,"Faulty","LowMedRD",paste0("LowRval(<",RvalTh,")")), values=c(poolCexList,4,3,5))+
          scale_colour_gradient("P-value", limits=c(0, 1), low="red", high="grey40")+
          geom_hline(yintercept=0,linetype="dashed", size=1)+
          geom_hline(yintercept=log2(delTh),linetype="dashed", size=0.5)+geom_hline(yintercept=log2(dupTh),linetype="dashed", size=0.5)+ 
	  geom_hline(yintercept=log2(1.5),linetype="solid", size=0.2)+geom_hline(yintercept=log2(0.5),linetype="solid", size=0.2)+
          annotate("text", x=dat$xvalues, y=MINV-0.2, label=dat$chrlabels, angle=90, size=2)
        
        ggsave(filename=paste0(plotdir,"/",MQ,"_dupdelTh",dupdelTh,"_",inSample,"_AllGene.pdf"), plot=last_plot(), width = 170, units= "mm") 
      }
    }
  }
  
  
  GenePlot=TRUE
  if(GenePlot==TRUE){
    ### For gene ##########
    geneList<-levels(RDTable$Gene)
    runTime<-round(NROW(geneList)/rangeX)  
    for(i in 0:(runTime-1)){
      geneList.run<-geneList[seq(from=i*rangeX+1, to=(i+1)*rangeX, by=1)]
      geneLen.run<-rangeX-length(geneList.run[geneList.run=="NA"])
      geneList.run<-geneList.run[1:geneLen.run]
      RDTable.run<-c()
      for(gene.run in geneList.run){
        RDTable.run<-rbind(RDTable.run,RDTable[RDTable$Gene==gene.run,])
      }      
      
      medCov<-unique(RDTable[RDTable$Type=="MedianRD",][c("Pool",inSample)])
      colnames(medCov)<-c("Pool","MedianRD")
      
      ampLen<-length(unique(RDTable.run$Amplicon_ID))
      
      xvalues<-rep(NA, ampLen)
      CNMvalues<-rep(NA, ampLen)
      CNUvalues<-rep(NA, ampLen)
      CNLvalues<-rep(NA, ampLen)
      pvalues<-rep(NA, ampLen)
      rvalues<-rep(NA, ampLen)
      cnvtypes<-rep(NA, ampLen)
      poolvalues<-rep(19, ampLen)
      exonvalues<-rep(19, ampLen)
      exonRvalues<-rep(19, ampLen)
      posvalues<-rep(19, ampLen)
      alphavalues<-rep(1.0, ampLen)
      
      labelX<-rep(NA,geneLen.run)
      labelTXT<-rep(NA,geneLen.run)
      
      xstart<-500
      start<-1
      for(x in 1:geneLen.run){
        gene=geneList.run[x]
        RDTable.gene<-RDTable.run[RDTable.run$Gene==gene,]
        ampLen.gene<-length(unique(RDTable.gene$Amplicon_ID))
        end<-start+ampLen.gene-1
        xend<-xstart+10*(ampLen.gene-1)
        
        xvalues[start:end]<-seq(from=xstart, to=xend, by=10)
        CNMvalues[start:end]<-as.numeric(as.matrix(RDTable.gene[RDTable.gene$Type=="CN_M",inSample]))
        CNUvalues[start:end]<-as.numeric(as.matrix(RDTable.gene[RDTable.gene$Type=="CI_U",inSample]))
        CNLvalues[start:end]<-as.numeric(as.matrix(RDTable.gene[RDTable.gene$Type=="CI_L",inSample]))
        pvalues[start:end]<-as.numeric(as.matrix(RDTable.gene[RDTable.gene$Type=="Pvalue",inSample]))
        rvalues[start:end]<-as.numeric(as.matrix(RDTable.gene[RDTable.gene$Type=="RegRvalue",inSample]))
        cnvtypes[start:end]<-as.character(as.matrix(RDTable.gene[RDTable.gene$Type=="CNVType",inSample]))
        poolvalues[start:end]<-as.character(as.matrix(RDTable.gene[RDTable.gene$Type=="CN_M","Pool"]))        
        posvalues[start:end]<-as.numeric(RDTable.gene[RDTable.gene$Type=="CN_M","Amplicon_Start"])
        
        exons<-RDTable.gene[RDTable.gene$Type=="CN_M","Exon"]
        exonvalues[start:end]<-as.character(as.matrix(gsub("Exon","E",exons)))
        exonRvalues[start:end]<-as.character(as.matrix(gsub("Exon","E",exons)))
        exon.now=""
        for(exonPos in c(start:end)){
          if(is.na(exonvalues[exonPos])|exonvalues[exonPos]==""){
            exonvalues[exonPos]<-""
          }else{
            exons<-unlist(strsplit(exonvalues[exonPos],"\\|"))
            if(as.character(exons[length(exons)])==exon.now){
              exonvalues[exonPos]<-""
            }else{
              out<-1
              for(i in c(1:length(exons))){
                if(as.character(exons[i])==exon.now){
                  out<-i+1
                }
              }
              exon<-paste(exons[c(out:length(exons))],collapse="|")
              exonvalues[exonPos]<-exon
              exon.now<-as.character(exons[length(exons)])
            }
          }
        }
        
        labelX[x]<-median(xvalues[start:end])
        TransList<-unique(as.character(RDTable.gene[RDTable.gene$Type=="CN_M","Transcript"]))
        TransID<-TransList[order(TransList)][length(TransList)]
        chr<-as.character(RDTable.gene[RDTable.gene$Type=="CN_M","Chr"])[1]
        if(TransID==""){
          labelTXT[x]<-paste0(gene,"\nchr",chr)
        }else{
          labelTXT[x]<-paste0(gene," (",as.character(TransID),")","\nchr",chr)
        }
        
        xstart<-xend+500
        start<-end+1
        x<-x+1
      }
      if(length(pvalues[is.na(pvalues)])==length(pvalues)){
        print("All p-values are nan!!! impossible to draw plot")
      }else{
        
        MAXV=3
        MINV=(-2)
        ###outlier handlin######
        temp.poolvalues<-poolvalues
        
        temp.CN_M<-CNMvalues
        temp.CN_M<-log2(temp.CN_M+0.000000001)
        temp.CN_M[temp.CN_M>MAXV]<-MAXV*1.0
        temp.CN_M[temp.CN_M<MINV]<-MINV*1.0
        CNMvalues<-temp.CN_M
        
        temp.CI_L<-CNLvalues
        temp.CI_L<-log2(temp.CI_L+0.000000001)
        temp.CI_L[temp.CI_L>MAXV]<-MAXV*1.0
        temp.CI_L[temp.CI_L<MINV]<-MINV*1.0
        CNLvalues<-temp.CI_L
        
        temp.CI_U<-CNUvalues
        temp.CI_U<-log2(temp.CI_U+0.000000001)  
        temp.CI_U[temp.CI_U>MAXV]<-MAXV*1.0
        temp.CI_U[temp.CI_U<MINV]<-MINV*1.0
        CNUvalues<-temp.CI_U 
       
        poolInfo<-c()
        refCexList<-c(15,16,17,18)
        poolCexList<-c()
        i<-1
        for(Pool in PoolList){
          poolCheck<-(medCov$Pool %in% c(Pool))
          if(length(poolCheck[poolCheck==TRUE])>0){
            if(as.numeric(as.matrix(medCov[medCov$Pool==Pool,2]))<50){
              temp.poolvalues[temp.poolvalues==Pool]<-"LowMedRD"
            }
          }
          poolInfo<-cbind(poolInfo, paste0(Pool,":",medCov[medCov$Pool==Pool,2]))
          poolCexList<-cbind(poolCexList,refCexList[i])
          i<-i+1
        }

        temp.poolvalues[rvalues<RvalTh]<-paste0("LowRval(<",RvalTh,")")
        temp.poolvalues[is.na(pvalues)]<-"Faulty"          

        poolvalues<-factor(temp.poolvalues,levels=c(PoolList,"Faulty","LowMedRD",paste0("LowRval(<",RvalTh,")")))    
        
        exonXvalues<-xvalues-5
        exonXvalues<-c(exonXvalues,exonXvalues[length(exonXvalues)]+10)
        exonList<-unique(unlist(strsplit(exonvalues,"\\|")))
        if(sum(exonList=="")!=length(exonList)){
          exonList<-exonList[exonList!=""]
          exonPosValues<-matrix(rep(NA, length(exonXvalues)*length(exonList)),length(exonList),length(exonXvalues))
          for(exon.now in exonList){
            check=TRUE
            for(i in c(1:length(exonList))){
              exonValue<-(-2.2)-(i*0.2)
              if(check==TRUE){
                if(min(grep(exon.now,exonRvalues))==1){
                  if(sum(is.na(exonPosValues[i,grep(exon.now,exonRvalues)]))==length(grep(exon.now,exonRvalues))){
                    pos<-grep(exon.now,exonRvalues)
                    exonPosValues[i,c(pos,pos[length(pos)]+1)]<-exonValue
                    check=FALSE
                  }
                }else if(sum(is.na(exonPosValues[i,c(min(grep(exon.now,exonRvalues))-1,grep(exon.now,exonRvalues))]))==length(grep(exon.now,exonRvalues))+1){
                  pos<-grep(exon.now,exonRvalues)
                  exonPosValues[i,c(pos,pos[length(pos)]+1)]<-exonValue
                  check=FALSE
                }
              }
            }
          }
          exonPosValues<-exonPosValues[rowSums(is.na(exonPosValues))!=length(exonPosValues[i,]), ]
        }else{
          exonPosValues<-t(rep(NA, length(exonXvalues)))
        }

	alphavalues[poolvalues=="Faulty"]<-0.2
        alphavalues[poolvalues=="LowMedRD"]<-0.6
        alphavalues[poolvalues=="LowRval"]<-0.6
        
        dat<-data.frame(xvalues, CNMvalues, CNUvalues,CNLvalues, pvalues, poolvalues, exonvalues, posvalues, alphavalues)
        dat2<-data.frame(exonXvalues,t(exonPosValues))
        
        medCovTXT=paste0(poolInfo,collapse=", ")
        titleTXT=paste0(inSample," (",MQ,")\n",medCovTXT)
        plot.now<-ggplot()+
          scale_x_continuous(breaks=c())+ 
          scale_y_continuous(limits=c(-3,3), breaks=seq(MINV,MAXV,1))+ 
          labs(x="Exons", y="Log2(observed read depth/expected read depth)")+theme_bw()+
          ggtitle(titleTXT)+
          theme(plot.title = element_text(size=20))+
          geom_point(data=dat, aes(x=dat$xvalues,y=dat$CNMvalues, shape=dat$poolvalues, color=dat$pvalues, alpha=dat$alphavalues), size=3)+
          geom_errorbar(data=dat, aes(x=dat$xvalues,ymin=dat$CNLvalues,ymax=dat$CNUvalues, color=dat$pvalues, alpha=dat$alphavalues))+
          scale_shape_manual("Pool",limits=c(PoolList,"Faulty","LowMedRD",paste0("LowRval(<",RvalTh,")")), values=c(poolCexList,4,3,5))+
          scale_colour_gradient("P-value", limits=c(0, 1), low="red", high="grey40")+
          scale_alpha(guide = 'none', limits=c(0,1))+ 
          geom_hline(yintercept=0,linetype="dashed", size=1)+
          geom_hline(yintercept=log2(delTh),linetype="dashed", size=0.5)+geom_hline(yintercept=log2(dupTh),linetype="dashed", size=0.5)+
	  geom_hline(yintercept=log2(1.5),linetype="solid", size=0.2)+geom_hline(yintercept=log2(0.5),linetype="solid", size=0.2)+
          annotate("text", x=labelX, y=MAXV-0.5, label=labelTXT, size=7)+
          annotate("text", x=xvalues, y=MINV-0.2, label=dat$exonvalues, angle=90, size=2)
        for(i in c(2:length(dat2))){
          if(sum(is.na(dat2[i]))!=NROW(dat2[i])){
            plot.now<-plot.now+geom_line(aes(x,y), data=data.frame(x=dat2$exonXvalues, y=unlist(dat2[i])), color="grey", size=3, alpha=0.5)
          }
        }    
        ggsave(filename=paste0(plotdir,"/",MQ,"_dupdelTh",dupdelTh,"_",inSample,"_",gene,".pdf"), plot=plot.now, width = 170, units= "mm")
      }
    }
  }
}
