args=commandArgs(T)

simulateSample<-function(depth,fetalConc,HetAltRatio,Perror,snpM,snpF,simulateCount,log){
  # Pat
  nP=0
  if(snpF>0){
     p0=fetalConc/2
  p1=Perror
  Pat=as.data.frame(x=t(replicate(simulateCount,simulatePath(snpF,depth,p0,p1))))
  colnames(Pat)=c("Ptrue","Pfalse","logPtrue","logPfalse")
  Pat$result=Pat$Ptrue>Pat$Pfalse
  write.table(Pat,file=paste0(log,".Pat.tsv"),sep="\t",row.names = F)
  #message(paste("Pat simulate finish",sum(Pat$result),sum(Pat$result)/simulateCount))
  nP=sum(Pat$result)
  }
 
  
  # Mat
  nM=0
  if(snpM>0){
  p0=HetAltRatio
  p1=HetAltRatio*(1-fetalConc)
  Mat=as.data.frame(x=t(replicate(simulateCount,simulatePath(snpM,depth,p0,p1))))
  colnames(Mat)=c("Ptrue","Pfalse","logPtrue","logPfalse")
  Mat$result=Mat$Ptrue>Mat$Pfalse
  write.table(Mat,file=paste0(log,".Mat.tsv"),sep="\t",row.names = F)
  nM=sum(Mat$result)
  }
  message(paste(log,depth,fetalConc,snpF,nP,snpM,nM))
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
HetAltRatio=0.48
Perror=0.02

count=as.numeric(args[1])

depth<-scan("workdir/input/depth.txt")
depth[depth==0]=200
ff<-scan("workdir/output/ff.txt")

result<-read.table("workdir/input/result.txt",stringsAsFactors = F,header = T)

info<-read.table("workdir/output/table1.txt",stringsAsFactors = F,header = F)
colnames(info)<-c("SID","Mp","Mn","Pp","Pn")
info$Depth=depth
info$Fra=ff
#info$snpM=info$Mp+info$Mn
#info$snpF=info$Pp+info$Pn

for(i in 1:nrow(info)){
	if(info[i,]$Mp*info[i,]$Mn>0){
		snpM=info[i,]$Mp*result$Mat[i]+info[i,]$Mn*(1-result$Mat[i])
	}else{
		snpM=info[i,]$Mp+info[i,]$Mn
	}
	if(info[i,]$Pp*info[i,]$Pn>0){
		snpF=info[i,]$Pp*result$Pat[i]+info[i,]$Pn*(1-result$Pat[i])
	}else{
		snpF=info[i,]$Pp+info[i,]$Pn
	}
	simulateSample(info[i,]$Depth,info[i,]$Fra,HetAltRatio,Perror,snpM,snpF,count,info[i,]$SID)
}
