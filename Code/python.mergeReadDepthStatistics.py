import sys

def getSampleList(inSampleInfoTxt):
  inFile=open(inSampleInfoTxt,'r')
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
  
  
def writeFile(sampleList, readDepthStatDir, outFileName):
  outFile=open(outFileName,'w')
  headCheck=False
  for sample in sampleList:
    sampleFileName=readDepthStatDir+"/"+sample+".readDepthStatistics.txt"
    inFile=open(sampleFileName,'r')
    inLineList=inFile.readlines()
    if(headCheck==False):
      headCheck=True
      outFile.write("".join(inLineList))
    else:
      outFile.write("".join(inLineList[1:]))
    inFile.close()
  outFile.close()
  
  
if __name__ == '__main__':
  inputs=list(sys.argv)
  outTag=inputs[1]
  inSampleInfoTxt=inputs[2]
  readDepthStatDir=inputs[3]
  
  sampleList=getSampleList(inSampleInfoTxt)
  outFileName=readDepthStatDir+"/"+outTag+".All.readDepthStatistics.txt"
  writeFile(sampleList, readDepthStatDir, outFileName)
  