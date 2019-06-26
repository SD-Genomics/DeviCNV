import sys
 
def getSampleList(inFileName):
  inFile=open(inFileName,'r')
  inLine=inFile.readline()
  headCheck=False
  sampleList=[]
  while(inLine):
    if(headCheck==False):
      headCheck=True
      headList=inLine.replace("\n","").replace("\r","").split("\t")
      sampleID=headList.index("Sample")
    else:
      inList=inLine.replace("\n","").replace("\r","").split("\t")
      sample=inList[sampleID]
      if(sample!="" and sample not in sampleList):
        sampleList.append(sample)
    inLine=inFile.readline()
  inFile.close()
  return sampleList
 
 
def mergeSampleFileList(batchTag, dupdelTh, MQ, CNVDir, sampleList):
  outFileName=CNVDir+"/"+batchTag+".CNV."+MQ+".dupdelTh"+dupdelTh+".txt"
  outFile=open(outFileName,'w')
  headCheck=False
  for sample in sampleList:
    qfileName=CNVDir+"/"+sample+".CNV."+MQ+".dupdelTh"+dupdelTh+".txt"
    try:
      qFile=open(qfileName,'r')
      qLine=qFile.readline()
      while(qLine):
        if(qLine.startswith("#")):
          if(headCheck==False):
            outFile.write(qLine)
        else:
          outFile.write(qLine)
        qLine=qFile.readline()
      headCheck=True
      qFile.close()   
    except:
     print("Unexpected error:", sys.exc_info()[0])
     print("Maybe all values are nan for "+sample)
         
  outFile.close()   
      
    
if __name__ == '__main__':
  inputs=list(sys.argv)
  batchTag=inputs[1]
  inSampleInfoTxt=inputs[2]
  CNVDir=inputs[3]
  MQList=inputs[4].split(",")
  dupdelList=inputs[5].split(",")
    
  sampleList=getSampleList(inSampleInfoTxt)
  for MQ in MQList:
    for dupdelTh in dupdelList:        
      mergeSampleFileList(batchTag, dupdelTh, MQ, CNVDir, sampleList)
