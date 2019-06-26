import numpy
import sys
from operator import itemgetter, attrgetter, methodcaller

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

    def __init__(self,ampID,chr,ampS,inS,inE,ampE,gene,trans,exon,pool):
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

    def allInfoList(self):
      return [self.ampID, self.chr, self.ampS,self.inS,self.inE,self.ampE, self.gene, self.trans, self.exon, self.pool]
        
    def head(self):
      return ["Amplicon_ID","Chr","Amplicon_Start","Insert_Start","Insert_End","Amplicon_End","Gene","Transcript","Exon","Pool"]


class RCPerSample:
    RCList=[]
    ampList=[]
    totalRC=0
    aveRC=0
    medRC=0
    stdRC=0
    sampleName=""
    gender=""
    MQ=""
    pool=""
    
    def __init__(self, RCList, ampList, sampleName, gender , MQ, pool):
        self.RCList=RCList[:]
        self.ampList=ampList     
        self.sampleName=sampleName
        self.gender=gender   
        self.MQ=MQ
        self.pool=pool
        self.totalRC=sum(self.RCList)
        self.aveRC=numpy.mean(self.RCList)
        self.medRC=numpy.median(self.RCList)
        self.stdRC=numpy.std(self.RCList)
        
    def normalizeGender(self):
        if(self.gender.lower()=="female"):
            oldtotalRC=self.totalRC
            for i in range(0, len(self.ampList)):
                if("X" in self.ampList[i].chr):
                    self.RCList[i]=self.RCList[i]/2  
            self.totalRC=sum(self.RCList) 
            self.aveRC=numpy.mean(self.RCList)
            self.medRC=numpy.median(self.RCList)      
            self.stdRC=numpy.std(self.RCList)     
            print(str(self.sampleName)+": Female normalization "+str(oldtotalRC)+" > "+str(self.totalRC))
        elif(self.gender.lower()=="male"):
          pass
        else:
          print(str(self.sampleName)+": wrongGender "+self.gender)
        
    def getRC(self, ampID):
        for pos in range(0, len(self.ampList)):
            amplicon=self.ampList[pos]
            if(amplicon.ampID==ampID):
                if(amplicon.chr=="Y" and self.gender=="female"):
                    return "Female"
                return self.RCList[pos]
        return "NA" 
    
    def showSampleInfo(self):
        return [self.sampleName, self.gender, self.MQ, self.pool]


def getSampleListAndDic(inSampleInfoTxt):
  inFile=open(inSampleInfoTxt,'r')
  inLine=inFile.readline()
  headCheck=False
  sampleList=[]
  sampleSexDic={}
  while(inLine):
    if(headCheck==False):
      headCheck=True
      headList=inLine.replace("\n","").replace("\r","").split("\t")
      sampleID=headList.index("Sample")
      sexID=headList.index("Sex")
    else:
      inList=inLine.replace("\n","").replace("\r","").split("\t")
      sample=inList[sampleID]
      sex=inList[sexID].lower()
      if(sample!="" and sample not in sampleList):
        sampleList.append(sample)
        sampleSexDic[sample]=sex
    inLine=inFile.readline()
  inFile.close()
  return [sampleList, sampleSexDic]


def getRCListPerSample(sampleList, readDepthDir, poolList, MQList, sampleSexDic):
  RCPerSampleDic={}
  ampliconDic={}
  for sample in sampleList:
    for pool in poolList:
      inFile=open(readDepthDir+sample+".readDepth."+pool+".txt")
      inLine=inFile.readline()
      headCheck=False
      RCDic={}
      ampList=[]
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
          ampliconDic[ampID]=Amplicon(ampID,chr,ampS,inS,inE,ampE,gene,exon,trans,pool)
          ampList.append(ampliconDic[ampID])
          for pos in range(0, len(inList)):
            if(header[pos].startswith("MQ") and header[pos] in MQList):
              MQ=header[pos]
              if(MQ not in RCDic):
                RCDic[MQ]=[]
              RCDic[MQ].append(float(inList[pos]))
        inLine=inFile.readline()
      inFile.close()
      for MQ in RCDic.keys():
        rcpersample=RCPerSample(RCDic[MQ], ampList, sample, sampleSexDic[sample], MQ, pool)
        rcpersample.normalizeGender()
        if(sample not in RCPerSampleDic):
          RCPerSampleDic[sample]=[]
        RCPerSampleDic[sample].append(rcpersample)
    
  return [RCPerSampleDic,ampliconDic]


def WriteNorRCFile(MQList, norRCFileName, sampleList, ampliconDic, RCPerSampleDic):
  head=ampliconDic[list(ampliconDic.keys())[0]].head()
  header=head+sampleList
  posID=header.index("Amplicon_Start")
  chrID=header.index("Chr")
  for MQ in MQList:
    outList=[]
    ampliconList=ampliconDic.keys()
    for ampliconID in ampliconList:
      amplicon=ampliconDic[ampliconID]
      outList.append(amplicon.allInfoList())
      for sampleName in sampleList:
        for rcpersample in RCPerSampleDic[sampleName]:
          if(rcpersample.MQ==MQ and rcpersample.pool==amplicon.pool):
            norRC=rcpersample.getRC(amplicon.ampID)
            outList[-1].append(norRC)
         
    outList=sorted(outList, key=itemgetter(posID))
    outList=sorted(outList, key=itemgetter(chrID))
    outFileName=norRCFileName+MQ+".txt"
    outFile=open(outFileName,'w')
    outFile.write("\t".join(numpy.array(header).astype(str))+"\n")
    for out in outList:
      outFile.write("\t".join(numpy.array(list(out)).astype(str))+"\n")
    outFile.close()


if __name__ == '__main__':
    inputs=list(sys.argv)
    batchTag=inputs[1]
    inSampleInfoTxt=inputs[2]
    readDepthDir=inputs[3]
    norRCDir=inputs[4]
    PoolList=inputs[5].split(",")
    MQList=inputs[6].split(",")

    [sampleList, sampleSexDic]=getSampleListAndDic(inSampleInfoTxt)
    [RCPerSampleDic,ampliconDic]=getRCListPerSample(sampleList, readDepthDir, PoolList, MQList, sampleSexDic)
    norRCFileName=norRCDir+"/"+batchTag+".readDepth.normalizedChrX."
    WriteNorRCFile(MQList, norRCFileName, sampleList, ampliconDic, RCPerSampleDic)

