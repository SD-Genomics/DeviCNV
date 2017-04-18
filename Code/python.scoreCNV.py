'''
Created on 2016. 8. 4.

@author: jino
'''
import sys
import numpy

def getHqSampleList(lqSampleFileName):
  inFile=open(lqSampleFileName,'r')
  inLine=inFile.readline()
  headCheck=False
  hqSampleList=[]
  sampleList=[]
  while(inLine):
    inList=inLine.replace("\n","").replace("\r","").split("\t")
    if(headCheck==False):
      headCheck=True
      sampleID=inList.index("Sample")
      lowQualID=inList.index("LowQualSample?")
    else:
      sample=inList[sampleID]
      lowQual=inList[lowQualID]
      if(lowQual=="HQ"):
        hqSampleList.append(sample)
      sampleList.append(sample)
    inLine=inFile.readline()
  inFile.close()
  return [sampleList, hqSampleList]
  

def getFilter(inFileName):
  inFile=open(inFileName,'r')
  inLine=inFile.readline()
  delfilterDic={}
  dupfilterDic={}
  headCheck=False
  while(inLine):
    if(headCheck==False):
      headCheck=True
      headList=inLine.lower().replace("\n","").replace("\r","").replace("#","").split("\t")
      attID=headList.index("attribute")
      delID=headList.index("delfilter")
      dupID=headList.index("dupfilter")          
    else:
      inList=inLine.lower().replace("\n","").replace("\r","").split("\t")
      if(inList[delID]!=""):
        if(inList[attID] not in delfilterDic):
          delfilterDic[inList[attID]]=[]
        delfilterDic[inList[attID]].append(inList[delID])
      if(inList[dupID]!=""):
        if(inList[attID] not in dupfilterDic):
          dupfilterDic[inList[attID]]=[]
        dupfilterDic[inList[attID]].append(inList[dupID])
    inLine=inFile.readline()
  inFile.close()
  filterList="#FilterList\t(Del):"
  filterAttList=list(delfilterDic.keys())
  filterAttList.sort()
  for filterAtt in filterAttList:
    filterList+=" "+filterAtt+str(delfilterDic[filterAtt])
  filterList+="\n#FilterList\t(Dup):"
  filterAttList=list(dupfilterDic.keys())
  filterAttList.sort()
  for filterAtt in filterAttList:
    filterList+=" "+filterAtt+str(dupfilterDic[filterAtt])
  return [delfilterDic, dupfilterDic, [filterList+"\n"]]
  
  
class ScoreStat:
  delfilterDic={}
  dupfilterDic={}
  sampleDic={}
  scoreID=[]
  sampleListHQ=[]
  sampleList=[]
  delMedian=[] 
  dupMedian=[] 
  allMedian=[] 
  delMean=[] 
  dupMean=[] 
  allMean=[] 
  delSum=[]
  dupSum=[]
  allSum=[]
  filterLen=0
  
  def __init__(self, delfilterDic, dupfilerDic):
    self.delfilterDic=delfilterDic
    self.dupfilterDic=dupfilterDic
    self.filterLen= len(delfilterDic)+2
    for i in range(0, len(delfilterDic)+1):
      self.scoreID.append(len(delfilterDic)-i)
    self.scoreID.append("raw")
    
  def putSampleListHQ (self, sampleListHQ):
    self.sampleListHQ =sampleListHQ
    
  def makeSampleScore (self, sampleList):
    self.sampleList=sampleList
    for sample in sampleList:
      self.sampleDic[sample]=[[0]*self.filterLen, [0]*self.filterLen, [0]*self.filterLen]
      
  def putSampleScore(self, sample, cnvType, score):
    if(cnvType=="del"):      
      self.sampleDic[sample][0][self.scoreID.index(score)]+=1
      self.sampleDic[sample][0][self.scoreID.index("raw")]+=1
    elif(cnvType=="dup"):
      self.sampleDic[sample][1][self.scoreID.index(score)]+=1
      self.sampleDic[sample][1][self.scoreID.index("raw")]+=1
    ## All #########################
    self.sampleDic[sample][2][self.scoreID.index(score)]+=1
    self.sampleDic[sample][2][self.scoreID.index("raw")]+=1
    
  def getMeanMedianScore(self):
    allScore=[[],[],[]]
    for i in range(self.filterLen):
      allScore[0].append([])
      allScore[1].append([])
      allScore[2].append([])
    for sample in self.sampleList:
      if(sample in self.sampleListHQ):
        for i in range(0,self.filterLen):
          allScore[0][i].append(self.sampleDic[sample][0][i])
          allScore[1][i].append(self.sampleDic[sample][1][i])
          allScore[2][i].append(self.sampleDic[sample][2][i])

    for allscores in allScore[0]:
      self.delMedian.append(numpy.median(allscores))
      self.delMean.append(numpy.around(numpy.mean(allscores),decimals=1))
      self.delSum.append(numpy.sum(allscores))
    for allscores in allScore[1]:
      self.dupMedian.append(numpy.median(allscores))   
      self.dupMean.append(numpy.around(numpy.mean(allscores),decimals=1))  
      self.dupSum.append(numpy.sum(allscores))   
    for allscores in allScore[2]:
      self.allMedian.append(numpy.median(allscores)) 
      self.allMean.append(numpy.around(numpy.mean(allscores),decimals=1)) 
      self.allSum.append(numpy.sum(allscores)) 

  def showStat(self):
    outAlls=[]
    outAlls.append(["##Score","(del):"]+self.scoreID+["(dup):"]+self.scoreID+["(total):"]+self.scoreID)
    outAlls.append(["##Sum(onlyHQSample)","(del):"]+self.delSum+["(dup):"]+self.dupSum+["(total):"]+self.allSum)
    outAlls.append(["##Mean(onlyHQSample)","(del):"]+self.delMean+["(dup):"]+self.dupMean+["(total):"]+self.allMean)
    outAlls.append(["##Median(onlyHQSample)","(del):"]+self.delMedian+["(dup):"]+self.dupMedian+["(total):"]+self.allMedian)
    outHQs=[]
    outLQs=[]
    for sample in self.sampleList:
      if(sample in self.sampleListHQ):
        outHQs.append(["#"+sample,"(del):"]+self.sampleDic[sample][0]+["(dup):"]+self.sampleDic[sample][1]+["(total):"]+self.sampleDic[sample][2])
      else:
        outLQs.append(["#(LQ)"+sample,"(del):"]+self.sampleDic[sample][0]+["(dup):"]+self.sampleDic[sample][1]+["(total):"]+self.sampleDic[sample][2])
    outList=[]
    for out in outAlls+outHQs+outLQs:
      outList.append("\t".join(numpy.array(out).astype(str))+"\n")
    return outList
        
        
def readCNVFile(inCNVFileName, delfilterDic, dupfilterDic, scoreStat):
  inFile=open(inCNVFileName,'r')
  inLine=inFile.readline()
  outHList=[]
  outVDic={}
  unifiedRegionIDList=[]
  regionIDDic={}
  headCheck=False
  while(inLine):
    if(inLine.startswith("#")):
      outHList.append(inLine)
      headLine=inLine
    else:
      if(headCheck==False):
        headCheck=True
        outHList[-1]="#Score\t"+outHList[-1].replace("#","")
        headList=headLine.lower().replace("\n","").replace("\r","").replace("#","").split("\t")
        regionIDID=headList.index("regionid")
        cnvTypeID=headList.index("cnvtype")
        sampleID=headList.index("sample")
        for headPos in xrange(len(headList)):
          if(headList[headPos] in delfilterDic):
            delfilterDic[headList[headPos]].append(headPos)
          if(headList[headPos] in dupfilterDic):
            dupfilterDic[headList[headPos]].append(headPos)    
        print(delfilterDic)
        print(dupfilterDic)
        print(headList)       
      ##########################################
      inList=inLine.replace("\n","").replace("\r","").split("\t")
      regionID=inList[regionIDID]
      ID=regionID.split(":")[-1]
      nameID=regionID.split(":")[0]+":"
      if("-sub" in ID):
        unifiedID=nameID+ID.split("-sub")[0]
        if(unifiedID not in regionIDDic):
          regionIDDic[unifiedID]=[]
        regionIDDic[unifiedID].append(regionID)
        unifiedRegionIDList.append(unifiedID)
      else:
        unifiedRegionIDList.append(regionID)
      sample=inList[sampleID] 
      cnvType=inList[cnvTypeID]  
      if(cnvType=="del"):
        filterDic=delfilterDic
      if(cnvType=="dup"):
        filterDic=dupfilterDic  
      passCnt=0
      for filterAtt in filterDic.keys():
        filterVList=filterDic[filterAtt][:-1]
        headPos=filterDic[filterAtt][-1]
        for filterV in filterVList:
          if(">=" in filterV):
            if(float(inList[headPos]) >= float(filterV.replace(">=",""))):
              passCnt+=1
          if((">" in filterV) and ("=" not in filterV)):
            if(float(inList[headPos]) > float(filterV.replace(">",""))):
              passCnt+=1
          if("<=" in filterV):
            if(float(inList[headPos]) <= float(filterV.replace("<=",""))):
              passCnt+=1
          if(("<" in filterV) and ("=" not in filterV)):
            if(float(inList[headPos]) < float(filterV.replace("<",""))):
              passCnt+=1     
      outVDic[regionID]=[passCnt, cnvType, sample, str(passCnt)+"\t"+inLine]
    inLine=inFile.readline()
  inFile.close()
  
  outVsList=[]
  unifiedRegionIDList=list(set(unifiedRegionIDList))
  for unifiedRegionID in unifiedRegionIDList:
    if(unifiedRegionID in regionIDDic):
      selectedVList=[]
      for regionID in regionIDDic[unifiedRegionID]:
        if(unifiedRegionID not in outVDic):
          selectedVList.append(outVDic[regionID])
        elif(outVDic[regionID][0]>outVDic[unifiedRegionID][0]):
          selectedVList.append(outVDic[regionID])
      outVsList+=selectedVList
      if(len(selectedVList)!=len(regionIDDic[unifiedRegionID])): 
        outVsList.append(outVDic[unifiedRegionID])
    else:
      outVsList.append(outVDic[unifiedRegionID])
  
  outVsList.sort(reverse=True)
  outVList=[]
  for outVs in outVsList:
    [score, cnvType, sample, results]=outVs
    outVList.append(results)
    scoreStat.putSampleScore(sample, cnvType, score)
  
  scoreStat.getMeanMedianScore()
  scoreStats=scoreStat.showStat()
      
  return [scoreStats,outHList,outVList]
  
  
def writeCNVFile(outCNVFileName, outList):
  outFile=open(outCNVFileName, 'w')
  for out in outList:
    outFile.write(out)
  outFile.close()
  
    
if __name__ == '__main__':
    inputs=list(sys.argv)
    batchTag=inputs[1]
    CNVDir=inputs[2]
    lqSampleDir=inputs[3]
    MQList=inputs[4].split(",")
    dupdelList=inputs[5].split(",")
    inScoringThTxt=inputs[6]

    for MQ in MQList:
      for dupdelTh in dupdelList:
        ### HighQualitySampleList ###########################
        lqSampleFileName_MQ=lqSampleDir+"/"+batchTag+".All.lowQualitySampleTest."+MQ+".txt"
        [sampleList, hqSampleList]=getHqSampleList(lqSampleFileName_MQ)
        ##########################################################
        inCNVFileName=CNVDir+"/"+batchTag+".All.CNV."+MQ+".dupdelTh"+dupdelTh+".txt"
        outCNVFileName=CNVDir+"/"+batchTag+".All.CNVWithScore."+MQ+".dupdelTh"+dupdelTh+".txt"
        [delfilterDic, dupfilterDic, filterList]=getFilter(inScoringThTxt)
        scoreStat=ScoreStat(delfilterDic, dupfilterDic)
        scoreStat.putSampleListHQ(hqSampleList)
        scoreStat.makeSampleScore (sampleList)
        [scoreStats,outHList,outVList]=readCNVFile(inCNVFileName, delfilterDic, dupfilterDic,scoreStat)
        writeCNVFile(outCNVFileName, scoreStats+filterList+outHList+outVList)
