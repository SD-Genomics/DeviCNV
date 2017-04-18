import numpy
from scipy import stats
import sys

def getHqSampleList(hqSampleFileName):
  inFile=open(hqSampleFileName,'r')
  inLine=inFile.readline()
  headCheck=False
  hqSampleList=[]
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
    inLine=inFile.readline()
  inFile.close()
  print("High-quality samples: "+str(hqSampleList))
  return hqSampleList


def getAmpRDDic(inFileName):
  inFile=open(inFileName,'r')
  inLine=inFile.readline()
  headCheck=False
  ampRDDic={}
  ampDic={}
  ampList=[]
  while(inLine):
    inList=inLine.replace("\n","").replace("\r","").split("\t")
    if(headCheck==False):
      headCheck=True
      ampliconIDID=inList.index("Amplicon_ID")
      poolID=inList.index("Pool")
      sampleID=poolID+1
      sampleList=inList[sampleID:] 
      headList=inList[:poolID+1]   
    else:
      ampliconID=tuple(inList[:poolID+1])
      pool=inList[poolID]
      ampDic[ampliconID]=pool
      sampleRDList=inList[sampleID:]
      ampRDDic[ampliconID]=sampleRDList
      ampList.append(ampliconID)
    inLine=inFile.readline()
  inFile.close()
  return [headList, sampleList, ampList, ampRDDic, ampDic]
  
  
def getMedCovDic(inFileName):
  inFile=open(inFileName,'r')
  inLine=inFile.readline()
  headCheck=False
  medCovDic={}
  medCovSampleList=[]
  while(inLine):
    inList=inLine.replace("\n","").replace("\r","").split("\t")
    if(headCheck==False):
      headCheck=True
      sampleID=inList.index("Sample")
      PoolID=inList.index("Pool")
      MQID=inList.index("MQ")
      medianID=inList.index("Median")      
    else:
      ID=(inList[PoolID], inList[MQID])
      if(ID not in medCovDic):
        medCovDic[ID]=[]
      medCovDic[ID].append(inList[medianID])
      if(inList[sampleID] not in medCovSampleList):
        medCovSampleList.append(inList[sampleID])
    inLine=inFile.readline()
  inFile.close()
  print("All samples")
  print(medCovSampleList)
  return [medCovSampleList, medCovDic]
    
  
def runBootstriping(sampleList, ampRDDic, medCovSampleList, medCovDic, ampliconID, Pool, MQ, dupTh, delTh, hqSampleList):
  faultyAmp=False
  ### filter out low quality sample ################################################
  x_raw=[]
  y_raw=[]  
  x_raw_highQ=[]
  y_raw_highQ=[]
  sampleList_highQ=[]
  for sampleName in sampleList:
    x_raw.append(float(medCovDic[(Pool, MQ)][medCovSampleList.index(sampleName)]))
    y_raw.append(float(ampRDDic[ampliconID][sampleList.index(sampleName)]))
    if(sampleName in hqSampleList):
      x_raw_highQ.append(float(medCovDic[(Pool, MQ)][medCovSampleList.index(sampleName)]))
      y_raw_highQ.append(float(ampRDDic[ampliconID][sampleList.index(sampleName)]))      
      sampleList_highQ.append(sampleName)
  ################## make regression model ################################################   
  slope, intercept, r_value, p_value, std_err = stats.linregress(x_raw_highQ,y_raw_highQ)
  ################## filter-out 20% outliner ##############################################
  y_raw_highQ_predict=intercept+slope*numpy.array(x_raw_highQ)
  residList=abs(numpy.array(y_raw_highQ)-y_raw_highQ_predict)
  Q3=numpy.percentile(residList,75)
  Q1=numpy.percentile(residList,25)
  IQR=Q3-Q1
  upperTh=Q3+1.5*IQR
  #maxThs=numpy.percentile(residList,90)
  filterOutPercent=0.2
  maxThs=residList[numpy.argsort(residList)][-1*round(len(residList)*filterOutPercent)]
  parsedSampleList=[]
  filterOutSampleList=[]
  for i in range(0, len(residList)):
    if((residList[i]<maxThs) or len(filterOutSampleList)>=round(len(residList)*filterOutPercent)):
      parsedSampleList.append(sampleList_highQ[i]) 
    elif(residList[i]>=maxThs):
      filterOutSampleList.append(sampleList_highQ[i])
  ########### print 2 type of warining #####################################################    
  if(len(sampleList_highQ)<10):
    print("Warning(Less 10 samples have >=50 medianCoverage)\t",ampliconID, Pool, MQ, "highQSampleList",str(sampleList_highQ),str(x_raw_highQ))
  if(len(filterOutSampleList)>round(len(residList)*filterOutPercent)):
    print("Warning(Over 10% samples have sample max residual)\t",ampliconID, Pool, MQ, maxThs,"filterOutSampleList",str(filterOutSampleList),str(residList))
  print(str(ampliconID)+": Total "+str(len(parsedSampleList))+" samples for bootstrapping")
  ########### bootstraping #################################################################  
  RUNTIME=1000
  yRatioList=[] 
  yPredictList=[]
  r_valueList=[]
  for t in range(0,RUNTIME):
    selectedSampleList=numpy.random.choice(parsedSampleList, len(parsedSampleList), replace=True)
    x=[]
    y=[]
    for sampleName in selectedSampleList:
      x.append(float(medCovDic[(Pool, MQ)][medCovSampleList.index(sampleName)]))
      y.append(float(ampRDDic[ampliconID][sampleList.index(sampleName)]))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    y_predict=intercept+numpy.array(x_raw)*slope
    yRatio=y_raw/y_predict
    yPredictList.append(list(y_predict))
    yRatioList.append(list(yRatio))
    r_valueList.append(r_value)
    ######## mark low quality amplicon #####
    if(slope<=0): 
      faultyAmp=True
  ########### mark low quality sample ############################################
  yPredictList=numpy.array(yPredictList)
  yRatioList=numpy.array(yRatioList)
  for i in range(0, len(yPredictList[0,:])):
      if(list(yPredictList[:,i]<0).count(True)>0):
          yPredictList[:,i]=["nan"]*len(yPredictList[:,i])
          yRatioList[:,i]=["nan"]*len(yRatioList[:,i])
  ############ write outList ###################################################
  x_raw=numpy.round(x_raw,3)
  y_raw=numpy.round(y_raw,3)
  yRatio_bottom=numpy.round(numpy.percentile(yRatioList,2.5, axis=0),3)
  yRatio_median=numpy.round(numpy.percentile(yRatioList,50,axis=0),3)
  yRatio_top=numpy.round(numpy.percentile(yRatioList,97.5, axis=0),3)
  yPredict_bottom=numpy.round(numpy.percentile(yPredictList,2.5, axis=0),3)
  yPredict_median=numpy.round(numpy.percentile(yPredictList,50,axis=0),3)
  yPredict_top=numpy.round(numpy.percentile(yPredictList,97.5, axis=0),3)
  ampP=[str(round(numpy.mean(r_valueList),2))]*len(x_raw) 
  yPredictList_t=numpy.transpose(yPredictList)
  yRatioList_t=numpy.transpose(yRatioList)
  X=[]
  Y=[]
  Y_M=[]
  Y_L=[]
  Y_U=[]
  CN_M=[]
  CI_L=[]
  CI_U=[]
  dupP=[]
  delP=[]
  P=[]
  CNVTYPE=[]
  for i in range(0,len(y_raw)): 
    ### calculate p-values for dup/del ########################################  
    dupPvalue=round(1-float(list(yRatioList_t[i]>dupTh).count(True))/RUNTIME,3)
    delPvalue=round(1-float(list(yRatioList_t[i]<delTh).count(True))/RUNTIME,3)  
    Pvalue=min(dupPvalue,delPvalue) 
    if((numpy.median(yPredictList_t[i])<y_raw[i]) and Pvalue<0.05):
        cnvType="dup"
    elif((numpy.median(yPredictList_t[i])>y_raw[i]) and Pvalue<0.05):
        cnvType="del"
    else:
        cnvType="neutral"
    ### mark low quality sample ################
    if(list(numpy.isnan(yRatioList_t[i])).count(True)>0):
        dupPvalue="nan"
        delPvalue="nan"
        Pvalue="nan"
        cnvType="faultySample"
    #### write out List #######################
    X.append(str(x_raw[i]))
    Y.append(str(y_raw[i]))
    Y_M.append(str(yPredict_median[i]))
    Y_L.append(str(yPredict_bottom[i]))
    Y_U.append(str(yPredict_top[i]))
    CN_M.append(str(yRatio_median[i]))
    CI_L.append(str(yRatio_bottom[i]))
    CI_U.append(str(yRatio_top[i]))
    dupP.append(str(dupPvalue))
    delP.append(str(delPvalue))
    P.append(str(Pvalue))
    CNVTYPE.append(cnvType)
  ### mark low quality amplicon #############  
  if(faultyAmp==True):
    dupP=["nan"]*len(dupP)
    delP=["nan"]*len(delP)
    P=["nan"]*len(P)
    CNVTYPE=["faultyAmp"]*len(CNVTYPE)
  return [X, Y, Y_M, Y_L, Y_U, CN_M, CI_L, CI_U, dupP, delP, P, CNVTYPE, ampP]
    
    
def writeAmpLinearRegression(headList, sampleList, ampList, ampRDDic, medCovSampleList, medCovDic, ampDic, modelFile, MQ, dupTh, delTh, hqSampleList):
  outFile=open(modelFile,'w')
  AmpliconHeaderList=headList+["MQ","Type"]
  outFile.write("\t".join(AmpliconHeaderList+sampleList)+"\n")
  for ampliconID in ampList:
    Pool=ampDic[ampliconID]
    try:
        [X, Y, Y_M, Y_L, Y_U, CN_M, CI_L, CI_U, dupP, delP, P, CNVTYPE,ampP]=runBootstriping(sampleList, ampRDDic, medCovSampleList, medCovDic, ampliconID, Pool, MQ, dupTh, delTh, hqSampleList)
        outFile.write("\t".join(list(ampliconID)+[MQ, "MedianRD"]+X)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "Y"]+Y)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "Y_L"]+Y_L)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "Y_M"]+Y_M)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "Y_U"]+Y_U)+"\n")        
        outFile.write("\t".join(list(ampliconID)+[MQ, "CI_L"]+CI_L)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "CN_M"]+CN_M)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "CI_U"]+CI_U)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "DupPvalue"]+dupP)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "DelPvalue"]+delP)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "Pvalue"]+P)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "CNVType"]+CNVTYPE)+"\n")
        outFile.write("\t".join(list(ampliconID)+[MQ, "RegRvalue"]+ampP)+"\n")
    except:
        print("error", ampliconID)
  outFile.close()


if __name__ == '__main__':
    inputs=list(sys.argv)    
    batchTag=inputs[1]
    readCountStatDir=inputs[2]
    norReadDepthDir=inputs[3]
    lqSampleDir=inputs[4]
    RDRatioDir=inputs[5]
    MQList=inputs[6].split(",")
    dupdelList=inputs[7].split(",")

    
    for MQ in MQList:
      hqSampleFileName_MQ=lqSampleDir+"/"+batchTag+".All.lowQualitySampleTest."+MQ+".txt"
      hqSampleList=getHqSampleList(hqSampleFileName_MQ)
      RDStaticFile=readCountStatDir+batchTag+".All.readDepthStatistics.txt"
      [medCovSampleList, medCovDic]=getMedCovDic(RDStaticFile)
      RDRatioFileName=norReadDepthDir+batchTag+".readDepth.normalizedChrX."+MQ+".txt"
      [headList, sampleList, ampList, ampRDDic, ampDic]=getAmpRDDic(RDRatioFileName)
      for dupdelTh in dupdelList:
        dupTh=float(dupdelTh.split("_")[0])
        delTh=float(dupdelTh.split("_")[1])
        modelFile=RDRatioDir+"/"+batchTag+".readDepthRatioFromLRModel."+MQ+".dupdelTh"+dupdelTh+".txt"
        writeAmpLinearRegression(headList, sampleList, ampList, ampRDDic, medCovSampleList, medCovDic, ampDic, modelFile, MQ, dupTh, delTh, hqSampleList)
        
    print("finish")
