
simulateSample<-function(depth,fetalConc,HetAltRatio,Perror,snpCount,simulateCount,log){
  # Pat
  p0=fetalConc/2
  p1=Perror
  Pat=as.data.frame(x=t(replicate(simulateCount,simulatePath(snpCount,depth,p0,p1))))
  colnames(Pat)=c("Ptrue","Pfalse","logPtrue","logPfalse")
  Pat$result=Pat$Ptrue>Pat$Pfalse
  write.table(Pat,file=paste0(log,".Pat.tsv"),sep="\t",row.names = F)
  message(paste("Pat simulate finish",sum(Pat$result),sum(Pat$result)/simulateCount))
  
  # Mat
  p0=HetAltRatio
  p1=HetAltRatio*(1-fetalConc)
  Mat=as.data.frame(x=t(replicate(simulateCount,simulatePath(snpCount,depth,p0,p1))))
  colnames(Mat)=c("Ptrue","Pfalse","logPtrue","logPfalse")
  Mat$result=Mat$Ptrue>Mat$Pfalse
  write.table(Mat,file=paste0(log,".Mat.tsv"),sep="\t",row.names = F)
  message(paste("Mat simulate finish",sum(Mat$result),sum(Mat$result)/simulateCount))
}

simulatePath<-function(snpCount,depth,p0,p1){
  path=sample(c(0,1),snpCount,replace=T)
  Pem0=c()
  n=length(path)
  for(i in 1:n){
    Pem0[i]=simulateSnp(path[i],depth,p0,p1)
  }
  Pem1=1-Pem0
  Ptrue=cumprod(Pem0)[n]
  Pfalse=cumprod(Pem1)[n]
  logPtrue=cumsum(log(Pem0))[n]
  logPfalse=cumsum(log(Pem1))[n]
  c(Ptrue,Pfalse,logPtrue,logPfalse)
}

simulateSnp<-function(hidenStatus,depth,p0,p1){
  RC=rpois(1,depth)
  if(hidenStatus==0){
    Pem0=simulateSnp0(RC,p0,p1)
  }else if(hidenStatus==1){
    Pem0=simulateSnp1(RC,p0,p1)
  }
  return(Pem0)
}

simulateSnp0<-function(RC,p0,p1){
  rc0=rbinom(1,RC,p0)
  P0=dbinom(rc0,RC,p0)
  P1=dbinom(rc0,RC,p1)
  Pem0=P0/(P0+P1)
  return(Pem0)
}

simulateSnp1<-function(RC,p0,p1){
  rc1=rbinom(1,RC,p1)
  P0=dbinom(rc1,RC,p0)
  P1=dbinom(rc1,RC,p1)
  Pem0=P1/(P0+P1)
  return(Pem0)
}


# Init
depth=200
fetalConc=0.1
HetAltRatio=0.48
Perror=0.02

# Model
snpCount=10
simulateSample(depth,fetalConc,HetAltRatio,Perror,snpCount)