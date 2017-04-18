'''
Created on 2016. 8. 4.

@author: jino
'''
import sys
import numpy
from operator import itemgetter, attrgetter
import vcf
from intervaltree import Interval, IntervalTree
from intervaltree_bio import GenomeIntervalTree
import math


class InfoPerAmplicon:
  ampID=""
  chr=""
  ampS=0
  inS=0
  inE=0
  ampE=0
  gene=""
  trans=""
  exon=""
  pval=-1
  cnm=-1
  ciLen=-1
  regRval=-1
  filter=""
    
  def __init__(self,ampID, chr, ampS, inS, inE, ampE, gene, trans, exon, pval, cnm, ciLen, regRval, filter):
    self.ampID=ampID
    self.chr=chr.replace("chr","")
    self.ampS=int(ampS) 
    self.inS=int(inS) 
    self.inE=int(inE)
    self.ampE=int(ampE)
    self.gene=gene
    self.trans=trans
    self.exon=exon.split("|")
    self.pval=float(pval)
    self.cnm=float(cnm)
    self.ciLen=float(ciLen)
    self.regRval=float(regRval)
    self.filter=filter
    
  def getInS(self):
    return self.inS


class GenePerSample:
    gene=""
    trans=""
    sample=""
    medCov=[]
    exonDupDic={}
    exonDelDic={}
    geneDupList=[]
    geneDelList=[]
    dupList=[]
    delList=[]
    pTh=0.05
    vcf=""
    coveredExonCnt=0
    totalAmpCnt=0
    qType=""
    dupTree=""
    delTree=""
    chrom=""
    start=-1
    end=-1
    minAmpCntInRegion=0
    dupTh=0
    delTh=0
    
    def __init__(self, gene, trans, sample, pTh, minAmpCntInRegion, dupTh, delTh):
        self.gene=gene
        self.sample=sample
        self.trans=trans
        self.medCov=[]
        self.exonDupDic={}
        self.exonDelDic={}
        self.geneDupList=[]
        self.geneDelList=[]   
        self.dupList=[]
        self.delList=[] 
        self.pTh=float(pTh)
        self.vcf=""
        self.coveredExonCnt=0
        self.totalExonCnt=0
        self.totalAmpCnt=0
        self.qType=""
        self.dupTree=GenomeIntervalTree()
        self.delTree=GenomeIntervalTree()
        self.chrom=""
        self.start=-1
        self.end=-1
        self.minAmpCntInRegion=minAmpCntInRegion
        self.dupTh=dupTh
        self.delTh=delTh
        
    def putTrans (self, trans):
        self.trans=trans
        
    def putMedCov(self,medCov):
        self.medCov.append(medCov)
      
    def putPvalue(self, ampID, chr, ampS, inS, inE, ampE, gene, trans, exon, type, pval, cnm, ciLen, regRval):
        self.chrom=chr
        if(self.start==-1):
          self.start=int(inS)
        if(self.end==-1):
          self.end=int(inE)
        self.start=min(self.start,int(inS))
        self.end=max(self.end,int(inE))
        if((pval!="nan") and (ciLen!="nan") and (cnm!="nan")): 
            filter="filter-out" 
            if(float(pval)<self.pTh):
                filter="filter-in"
            info=InfoPerAmplicon(ampID, chr, ampS, inS, inE, ampE, gene, trans, exon, pval, cnm, ciLen, regRval, filter)
            if(type=="dup"):
                self.dupList.append(info)
                self.dupTree.addi(chr,int(inS),int(inS)+1,info)
            elif(type=="del"):
                self.delList.append(info)
                self.delTree.addi(chr,int(inS),int(inS)+1,info)

    def getCandidateRegion(self, qType, CBSTree):
        if(qType=="dup"):
            sortedList=sorted(self.dupList,key=attrgetter('inS'))
        elif(qType=="del"):
            sortedList=sorted(self.delList,key=attrgetter('inS'))
        ### correct Exon Info ##################   
        totalExonList=[] 
        for Info in sortedList:
            for exon in Info.exon:
              if("Exon" in exon):
                totalExonList.append(int(exon.replace("Exon","")))            
        totalExonList=list(set(totalExonList))
        self.coveredExonCnt=len(totalExonList)   
        self.totalAmpCnt=len(sortedList)  
        self.qType=qType
        ##################################################
        #### UnifiedRegions ##############################
        UnifiedRegion={}
        ID=1
        segmentList=CBSTree[self.chrom].search(self.start, self.end+1) 
        segmentList=sorted(segmentList,key=itemgetter(1))
        for segment in segmentList:
          [segChr, segStart, segEnd, segNo, segMean]=segment[-1]
          if(self.qType=="dup"):
            dupList=[]
            for overlapAmps in list(self.dupTree[segChr].search(segStart, segEnd)):
              dupList.append(overlapAmps[-1])
            if(dupList!=[] and len(dupList)>=minAmpCntInRegion):
              dupList=sorted(dupList,key=attrgetter('inS'))
              UnifiedRegion[self.sample+"_"+self.gene+"_"+"dup"+":"+str(ID)]=dupList  
          if(self.qType=="del"):
            delList=[]
            for overlapAmps in list(self.delTree[segChr].search(segStart, segEnd)):
              delList.append(overlapAmps[-1])
            if(delList!=[] and len(delList)>=minAmpCntInRegion):
              delList=sorted(delList,key=attrgetter('inS'))
              UnifiedRegion[self.sample+"_"+self.gene+"_"+"del"+":"+str(ID)]=delList
          ID+=1
        #### Regions in a unifiedRegions ##################
        UnifiedRegionIDList=list(UnifiedRegion.keys())
        for ID in UnifiedRegionIDList:
          subID=1 
          filterIn=[]      
          for Info in UnifiedRegion[ID]:
            if(Info.filter=="filter-in"):
              filterIn.append(Info)
            else:
              if(filterIn!=[]):
                if(len(filterIn)>=minAmpCntInRegion and filterIn!=UnifiedRegion[ID]):
                  UnifiedRegion[ID+"-sub"+str(subID)]=filterIn
                  subID+=1
              filterIn=[]
          if(filterIn!=[]):
            if(len(filterIn)>=minAmpCntInRegion and filterIn!=UnifiedRegion[ID]):
              UnifiedRegion[ID+"-sub"+str(subID)]=filterIn
        ###################################################    
        TotalUnifiedRegionIDList=list(UnifiedRegion.keys())
        TotalUnifiedRegionIDList.sort()
        regionList=[]
        for ID in TotalUnifiedRegionIDList:
          infos=UnifiedRegion[ID]
          infosAnnotation=self.getRegionInfo(infos)
          if(infosAnnotation!=""):
            regionList.append(ID+"\t"+infosAnnotation)
        return regionList
       
    def getRegionInfo(self, infos):
        ampCnt=len(infos)
        exons=[]
        pvals=[]
        amps=[]
        cnms=[]
        ciLens=[]
        regRvals=[]
        chr=infos[0].chr
        start=infos[0].inS
        end=infos[-1].inE
        filterCheck=[]
        for info in infos:
          amps+=[info.ampID]
          exons+=info.exon
          pvals+=[info.pval]
          cnms+=[info.cnm]
          ciLens+=[info.ciLen]
          regRvals+=[info.regRval]
          filterCheck+=[info.filter]
        for i in range(len(exons)):
          if(exons[i]==""):
            exons[i]="Intron"
        selectedExons=list(set(exons))
        if("Intron" in selectedExons):
          selectedExons.remove("Intron")
        selectedExonCnt=len(selectedExons)
        try:
          exonInRegionRatio=round(float(selectedExonCnt)/self.coveredExonCnt,2)
        except:
          exonInRegionRatio="nan"
        medCovList=list(set(self.medCov))
        medCovList.sort()      
        region=[self.qType, self.sample, ",".join(medCovList), chr, start, end, self.gene, self.trans]
        region+=[selectedExonCnt, self.coveredExonCnt, exonInRegionRatio , "|".join(exons)]
        region+=[ampCnt, self.totalAmpCnt ,"|".join(amps), ",".join(numpy.array(pvals).astype(str)),round(float(filterCheck.count("filter-in"))/len(filterCheck),2)]
        region+=[numpy.mean(cnms),numpy.std(cnms),",".join(numpy.array(cnms).astype(str)), numpy.mean(ciLens), numpy.mean(regRvals)]     
        
        if(self.qType=="dup" and numpy.mean(cnms)<=self.dupTh):
          return ""
        elif(self.qType=="del" and numpy.mean(cnms)>=self.delTh):
          return ""
        else:
          return "\t".join(numpy.array(region).astype(str))+"\n"
          
    def getHead(self):
      HeadList=["RegionID","CnvType","Sample","MedianRDOfSample","Chr","Start","End","Gene","Transcript","ExonCntInRegion","CoveredExonCnt","ExonInRegionRatio","Exons","AmpCntInRegion","TotalAmpCnt","Amplicons"]
      HeadList+=["Pvalues","FilterInAmpRatio","AverageOfReadDepthRatios","STDOfReadDepthRatios","ReadDepthRatios","AverageOfCIs","AverageOfR2vals"]
      return HeadList
      
       
def getRDRatioDic(inFileName, batchTag, pTh, sample, CBSTree, minAmpCntInRegion, dulTh, delTh):
    inFile=open(inFileName,'r')
    inLine=inFile.readline()
    headCheck=False
    genePerSampleDic={}
    medCovList=""
    while(inLine):
      if(headCheck==False):
        headCheck=True
        header=inLine.replace("\n","").split("\t")
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
        typeID=header.index("Type")
        sampleID=header.index(sample)
      else:        
        inList=inLine.replace("\n","").split("\t")
        type=inList[typeID]
        if(type=="MedianRD"):
            if(medCovList!=""):
              ampID=medCovList[ampIDID]
              chr=medCovList[chrID]
              ampS=medCovList[ampSID]
              inS=medCovList[inSID]
              inE=medCovList[inEID]
              ampE=medCovList[ampEID]
              gene=medCovList[geneID]
              trans=medCovList[transID]
              exon=medCovList[exonID]
              pool=medCovList[poolID]
              regRval=regRvalList[sampleID]
              medCov=pool+":"+medCovList[sampleID]
              cnm=cnmList[sampleID]
              try:
                ciLen=float(ciuList[sampleID])-float(cilList[sampleID])
              except:
                ciLen="nan"
              dupPval=dupPvalList[sampleID]
              delPval=delPvalList[sampleID]
              cnvtype=cnvtypeList[sampleID]
              if("low" in cnvtype):
                dupPval="nan"
                delPval="nan"
              ID=(gene,sample)
              if(ID not in genePerSampleDic):
                genePerSampleDic[ID]=GenePerSample(gene,trans, sample, pTh,minAmpCntInRegion, dulTh, delTh)
              if(trans!=""):
                genePerSampleDic[ID].putTrans(trans)
              genePerSampleDic[ID].putMedCov(medCov)
              genePerSampleDic[ID].putPvalue(ampID, chr, ampS, inS, inE, ampE, gene, trans, exon, "del", delPval, cnm, ciLen, regRval) 
              genePerSampleDic[ID].putPvalue(ampID, chr, ampS, inS, inE, ampE, gene, trans, exon, "dup", dupPval, cnm, ciLen, regRval)  
            else:
              pass #frist line          
            medCovList=inList
        elif(type=="CI_U"):
            ciuList=inList
        elif(type=="CI_L"):
            cilList=inList
        elif(type=="CN_M"):
            cnmList=inList
        elif(type=="DupPvalue"):
            dupPvalList=inList
        elif(type=="DelPvalue"):            
            delPvalList=inList
        elif(type=="CNVType"):
            cnvtypeList=inList
        elif(type=="RegRvalue"):
            regRvalList=inList
        inLine=inFile.readline()
    inFile.close()
    return genePerSampleDic


def getRegionCandidate(genePerSampleDic, outFileName, pTh, CBSTree,minAmpCntInRegion, dupTh, delTh, MQ):
    outFile=open(outFileName,'w')
    outFile.write("#P-value<"+str(pTh)+"\n")
    outFile.write("#Minimum-number of amplicons to extract small CNVs >="+str(minAmpCntInRegion)+"\n")
    outFile.write("#dupTh:"+str(dupTh)+", delTh:"+str(delTh)+"\n")
    outFile.write("#MQV>="+str(MQ)+"\n")
    outFile.write("#"+"\t".join(genePerSampleDic[list(genePerSampleDic.keys())[0]].getHead())+"\n")
    genePerSampleIDs=list(genePerSampleDic.keys())
    for genepersampleID in genePerSampleIDs:
        genepersample=genePerSampleDic[genepersampleID]
        out=genepersample.getCandidateRegion("del", CBSTree)
        out+=genepersample.getCandidateRegion("dup", CBSTree)
        if(len(out)>0):
            outFile.write("".join(out))
    outFile.close()
       
       
def readCBSFile(inFileName):
  inFile=open(inFileName,'r')
  inLine=inFile.readline()
  headCheck=False
  dupList=[]
  delList=[]
  CBSTree=GenomeIntervalTree()
  while(inLine):
    if(inLine.startswith("#")):
      pass
    else:
      if(headCheck==False):
        headCheck=True
        headList=inLine.replace("\n","").split("\t")
        chrID=headList.index("chromosome")
        startID=headList.index("start")
        endID=headList.index("end")
        noID=headList.index("nbrOfLoci")
        meanID=headList.index("mean")
      else:
        inList=inLine.replace("\n","").split("\t")
        chr=inList[chrID].replace("23","X").replace("24","Y").replace("25","MT")
        no=int(inList[noID])
        try:
          mean=float(inList[meanID])
          start=int(inList[startID])
          end=int(inList[endID])
        except:
          mean="nan"
          start=-1
          end=-1
        if(mean!="nan"):
          if(mean>0):
            dupList.append([inList[chrID], inList[startID], inList[endID], mean])
          if(mean<0):
            delList.append([inList[chrID], inList[startID], inList[endID], mean]) 
          CBSTree.addi(chr, start, end+1, [chr, start, end+1, no, mean])
    inLine=inFile.readline()
  inFile.close()
  return [dupList, delList, CBSTree]
    

if __name__ == '__main__':
  inputs=list(sys.argv)
  batchTag=inputs[1]
  inSample=inputs[2]
  RDRatioDir=inputs[3]
  CBSDir=inputs[4]
  CNVDir=inputs[5]
  MQList=inputs[6].split(",")
  dupdelList=inputs[7].split(",")
  
    
  pTh=0.05 ## p-value threshold
  minAmpCntInRegion=2 ## The minimum number of significant amplicons to extract small region CNV candidates
    
  for MQ in MQList:
    for dupdelTh in dupdelList:
      dupTh=float(dupdelTh.split("_")[0])
      delTh=float(dupdelTh.split("_")[1])
      inFileName=CBSDir+"/"+inSample+".CBS."+MQ+".dupdelTh."+dupdelTh+".tsv"
      outFileName=CNVDir+"/"+inSample+".CNV."+MQ+".dupdelTh"+dupdelTh+".txt"
      try:
        [dupList,delList,CBSTree]=readCBSFile(inFileName)
        inFileName2=RDRatioDir+batchTag+".readDepthRatioFromLRModel."+MQ+".dupdelTh"+dupdelTh+".txt"
        genePerSampleDic=getRDRatioDic(inFileName2, batchTag, pTh, inSample, CBSTree, minAmpCntInRegion, dupTh, delTh)
        getRegionCandidate(genePerSampleDic, outFileName, pTh, CBSTree,minAmpCntInRegion, dupTh, delTh, MQ)
      except:
          print("Unexpected error:", sys.exc_info()[0])
          print("Maybe all values are nan for "+inSample)