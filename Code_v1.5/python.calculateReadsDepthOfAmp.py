'''
Created on 2016. 5. 10.

@author: jino
'''
import random
import numpy
import pysam
import sys
from intervaltree import Interval, IntervalTree
from intervaltree_bio import GenomeIntervalTree

class Amplicon:
    ampID=""
    chr=""
    ampS=0
    inS=0
    inE=0
    ampE=0
    gene=""
    trans=""
    exon=""
    pool=""
    datType=""
    mappedReadList=[]
    readDepthList=[]
    
    def __init__(self,ampID,chr,ampS,inS,inE,ampE,gene,trans,exon,pool,datType):
        self.ampID=ampID
        self.chr=chr.replace("chr","")
        self.ampS=int(ampS) 
        self.inS=int(inS) 
        self.inE=int(inE)
        self.ampE=int(ampE)
        self.gene=gene
        self.trans=trans
        self.exon=exon
        self.pool=pool.split("_")[-1]
        self.datType=datType
        self.mappedReadList=[]
        self.readDepthList=[]
        
    def putMappedRead(self, inMappedRead):
        self.mappedReadList.append(inMappedRead)
    
    def getOverlapRatio(self, rChr, rS, rE):  
      if(rChr.replace("chr","")!=self.chr):
        return 0.0
      else:
        rLen=rE-rS
        overlapLen=min(rE, self.ampE)-max(rS, self.ampS)
        overlapRatio=float(overlapLen)/float(rLen)
        if(overlapRatio>1):
          return 1.0
        else:
          return overlapRatio 

    def getReadDepthPCR(self, MQList): 
        ampliconLength=self.inE-self.inS
        depthPerSiteDic={}
        for MQ in MQList:
            depthPerSiteDic[MQ]=[0]*ampliconLength                      
        for pos in range(0, ampliconLength):
            nowS=self.inS+pos
            for read in self.mappedReadList:
                if(read.pos<=nowS and nowS+1<=read.pos+read.alen):
                    for MQ in MQList:
                        if(read.mapq>=MQ):
                            depthPerSiteDic[MQ][pos]+=1                            
        readDepthOutList=[]
        for MQ in MQList:
          readDepth=0
          for read in self.mappedReadList:
            if(read.mapq>=MQ):
              readDepth+=1
          readDepthOutList.append(readDepth)   
        self.readDepthList=readDepthOutList
    
    def getReadDepthHYB(self, MQList):  ## using insert
        ampliconLength=self.inE-self.inS
        depthPerSiteDic={}
        for MQ in MQList:
            depthPerSiteDic[MQ]=[0]*ampliconLength                      
        for pos in range(0, ampliconLength):
            nowS=self.inS+pos
            for read in self.mappedReadList:
                if(read.pos<=nowS and nowS+1<=read.pos+read.alen):
                    for MQ in MQList:
                        if(read.mapq>=MQ):
                            depthPerSiteDic[MQ][pos]+=1                            
        readDepthOutList=[]
        for MQ in MQList:
          depCov=0
          for depth in depthPerSiteDic[MQ]:
            depCov+=depth
          readDepthOutList.append(round(float(depCov)/ampliconLength,3))                           
        self.readDepthList=readDepthOutList
                      
    def runGetReadDepth(self, MQList):
      if(self.datType=="HYB"):
        self.getReadDepthHYB(MQList)
      elif(self.datType=="PCR"):
        self.getReadDepthPCR(MQList)
      else:
        print(self.datType, "unknown data")
           
    def allInfoList(self):
      return [self.ampID, self.chr, self.ampS,self.inS,self.inE,self.ampE, self.gene, self.trans, self.exon, self.pool]
        
    def head(self):
      return ["Amplicon_ID","Chr","Amplicon_Start","Insert_Start","Insert_End","Amplicon_End","Gene","Transcript","Exon","Pool"]
        
        
def MakeAmpliconDic(inAmpliconTxt, datType):
    inFile=open(inAmpliconTxt, 'r')
    inLine=inFile.readline()
    ampliconDic={}
    ampLocalDic={}
    ampliconList=[]
    ampTree=GenomeIntervalTree()
    headCheck=False
    while(inLine):
      if(headCheck==False):
        headCheck=True
        header=inLine.replace("\n","").replace("\r","").split("\t")
        print(header)
        ampIDID=header.index("Amplicon_ID")
        chrID=header.index("Chr")
        ampSID=header.index("Amplicon_Start")
        inSID=header.index("Insert_Start")
        inEID=header.index("Insert_End")
        ampEID=header.index("Amplicon_End")
        geneID=header.index("Gene")
        transID=header.index("Transcript")
        exonID=header.index("Exon")
        poolID=header.index("Pool")     
      else:
        inList=inLine.replace("\n","").replace("\r","").split("\t")
        ampID=inList[ampIDID]
        chr=inList[chrID].replace("chr","")
        ampS=inList[ampSID]
        inS=int(inList[inSID])
        inE=int(inList[inEID])
        ampE=inList[ampEID]
        gene=inList[geneID]
        exon=inList[exonID]
        trans=inList[transID]
        pool=inList[poolID]
        if(ampID not in ampLocalDic):
            ampliconList.append(ampID)
            ampLocalDic[ampID]=Amplicon(ampID,chr,ampS,inS,inE,ampE,gene,exon,trans,pool,datType)
            ampTree.addi(chr,inS+1,inE+1,ampID) ## [start, end)
        else:
            print("Error!! : Not unique Amplicon_ID : "+ampID)
            break
      inLine=inFile.readline()
    inFile.close()    
    
    for ampliconID in ampliconList:
        amplicon=ampLocalDic[ampliconID]
        pool=amplicon.pool
        if(pool not in ampliconDic):
            ampliconDic[pool]=[]
        ampliconDic[pool].append(amplicon)
    print("Total Amplicons: "+str(len(ampLocalDic.keys())))
    print("ampTree made!")
    return [ampliconDic, ampTree]
   

def MapReadinBamPCR(inBamFile, ampliconDic, ampTree, dedupOp, MQList):
    ampliconList=[]
    poolList=list(ampliconDic.keys())
    poolList.sort()
    for pool in poolList:
        ampliconList+=ampliconDic[pool]
        
    inBam=pysam.Samfile(inBamFile,'rb')
    for read in inBam:
        if(read.is_unmapped):
            pass
        else:
            if(read.is_duplicate):
                if(dedupOp=="true"):
                    continue
            overlapAmpTreeList=ampTree[inBam.getrname(read.rname).replace("chr","")].search(read.pos+1, read.pos+read.alen+1) ## [start, end)
            if(len(overlapAmpTreeList)==0):
                pass
            else:
                overlapAmpIDList=[] 
                for overlapAmpTree in overlapAmpTreeList:
                  overlapAmpIDList.append(overlapAmpTree[-1])
                
                overlapAmpList=[]   
                for amplicon in ampliconList:
                  if(amplicon.ampID in overlapAmpIDList):
                    overlapAmpList.append(amplicon)

                overlapRatioList=[]   
                ampLenList=[]         
                for amplicon in overlapAmpList:
                  overlapRatioList.append(amplicon.getOverlapRatio(inBam.getrname(read.rname).replace("chr",""), read.pos, read.pos+read.alen)) 
                  ampLenList.append(amplicon.ampE-amplicon.ampS)

                maxValue=max(overlapRatioList)
                overlapAmpList2=[]
                overlapRatioList2=[]   
                ampLenList2=[]  
                for i in range(0,len(overlapAmpList)):
                    if(maxValue==overlapRatioList[i]):
                        overlapAmpList2.append(overlapAmpList[i])
                        overlapRatioList2.append(overlapRatioList[i])
                        ampLenList2.append(ampLenList[i])
                        
                minAmpLen=min(ampLenList2) 
                overlapAmpList3=[]
                overlapRatioList3=[]   
                ampLenList3=[]                              
                for j in range(0,len(overlapAmpList2)):  
                    if(minAmpLen==ampLenList2[j]):
                        overlapAmpList3.append(overlapAmpList2[j])
                        overlapRatioList3.append(overlapRatioList2[j])
                        ampLenList3.append(ampLenList2[j])          
                                  
                mappedAmp=overlapAmpList3[int((random.random()*10000))%(len(overlapAmpList3))]
                mappedAmp.mappedReadList.append(read)                   
    
    for amplicon in ampliconList:
      amplicon.runGetReadDepth(MQList)
                  
    return ampliconDic


def MapReadinBamHYB(inBamFile, ampliconDic, ampTree, dedupOp, MQList):
    ampliconList=[]
    poolList=list(ampliconDic.keys())
    poolList.sort()
    for pool in poolList:
        ampliconList+=ampliconDic[pool]
        print(pool)
        
    inBam=pysam.Samfile(inBamFile,'rb')
    
    for read in inBam:
        if(read.is_unmapped):
            pass
        else:
            if(read.is_duplicate):
                if(dedupOp=="true"):
                    continue
            overlapAmpTreeList=ampTree[inBam.getrname(read.rname).replace("chr","")].search(read.pos+1, read.pos+read.alen+1) ## [start, end)
            if(len(overlapAmpTreeList)==0):
                pass
            else:          
                overlapAmpIDList=[]
                for overlapAmpTree in overlapAmpTreeList:
                  overlapAmpIDList.append(overlapAmpTree[-1])
                for amplicon in ampliconList:
                  if(amplicon.ampID in overlapAmpIDList):
                    amplicon.mappedReadList.append(read)
                 
    for amplicon in ampliconList:
      amplicon.runGetReadDepth(MQList)
                  
    return ampliconDic
    
    
def WriteReadDepthFile(ampliconDic, outFileName, MQList):
    ### write file per pool ###########################
    ampliconList=list(ampliconDic.keys())
    ampliconList.sort()
    for pool in ampliconList:
        #### write attributes ##########################
        outFile=open(outFileName+"."+pool+".txt",'w')
        header=ampliconDic[pool][0].head()
        outFile.write("\t".join(header))
        for MQ in MQList:
            outFile.write("\tMQ"+str(MQ))
        outFile.write("\n")
        #### write values per amplicon ################
        for amplicon in ampliconDic[pool]:
            outFile.write("\t".join(numpy.array(amplicon.allInfoList()).astype(str)))
            readDepthOutList=amplicon.readDepthList
            outFile.write("\t"+"\t".join(numpy.array(readDepthOutList).astype(str)))
            outFile.write("\n")            
        outFile.close()


def WriteMappedReadDepthStatFile(ampliconDic, RCstaticFileName, MQList, inSample):
    staticFile=open(RCstaticFileName+".txt",'w')
    staticFile.write("Sample\tPool\tMQ\tMean\tMedian\tStandardDeviation\tSum\n")
    ### write file per pool ###########################
    ampliconList=list(ampliconDic.keys())
    ampliconList.sort()
    for pool in ampliconList:
        totalReadDepthOutList=[]
        for amplicon in ampliconDic[pool]:
            readDepthOutList=amplicon.readDepthList
            totalReadDepthOutList.append(readDepthOutList)  
        #### write StaticFile per Pool+MQ #############
        totalReadDepthOutList=numpy.transpose(totalReadDepthOutList)
        for i in range(0,len(MQList)):
            MQ=MQList[i]
            RCList=totalReadDepthOutList[i]
            staticList=[round(numpy.mean(RCList),2),round(numpy.median(RCList),2), round(numpy.std(RCList),2), round(numpy.sum(RCList))]
            staticFile.write(inSample+"\t"+pool+"\tMQ"+str(MQ)+"\t"+"\t".join(numpy.array(staticList).astype(str))+"\n")        
    #####################################################    
    staticFile.close()    
    
    
if __name__ == '__main__':
    inputs=list(sys.argv)
    inSample=inputs[1]
    inBamDir=inputs[2]
    inAmpliconTxt=inputs[3]
    readDepthDir=inputs[4]
    readDepthStatDir=inputs[5]
    dedupOp=inputs[6].lower()
    datType=inputs[7]
    MQList=list(numpy.array(inputs[8].replace("MQ","").split(",")).astype(int))

    [ampliconDic, ampTree]=MakeAmpliconDic(inAmpliconTxt,datType)
    
    inBamFile=inBamDir+inSample+".bam"     
    if(datType=="HYB"):
      ampliconDic=MapReadinBamHYB(inBamFile, ampliconDic, ampTree, dedupOp, MQList)
    elif(datType=="PCR"):
      ampliconDic=MapReadinBamPCR(inBamFile, ampliconDic, ampTree, dedupOp, MQList)
    else:
      print("ERROR !! Unknown data type")
      
    readDepthFile=readDepthDir+inSample+".readDepth"
    WriteReadDepthFile(ampliconDic, readDepthFile, MQList)
    RCStaticFile=readDepthStatDir+inSample+".readDepthStatistics"
    WriteMappedReadDepthStatFile(ampliconDic, RCStaticFile, MQList, inSample)



