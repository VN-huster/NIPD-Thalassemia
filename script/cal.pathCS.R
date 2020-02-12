#args=commandArgs(T)

sList=paste0("F",c(paste0(0,1:9),10:59))
fList=paste0(sList,"_nipt.fout")
mList=paste0(sList,"_nipt.mout")
info<-cbind(sList,fList,mList)
info<-data.frame(cbind(sList,fList,mList),stringsAsFactors=F)
colnames(info)<-c("SID","Ffile","Mfile")
info$mCount=0
info$mCS=0
info$fCount=0
info$fCS=0

#head(info)

for(i in 1:nrow(info)){
  message(info[i,]$SID)
  if(file.exists(file=info[i,]$Ffile)){
    read.table(info[i,]$Ffile,header=T,stringsAsFactors=F)->f
    info$fCount[i]=nrow(f)
    ptrue=f$P0
    pfalse=f$P1
    for(j in 1:nrow(f)){
	if((!is.null(f$fetal_F[j])&&f$fetal_F[j]==2)||(is.null(f$fetal_F[j])&&f$P0<0.5)){
	    ptrue[j]=f$P1[j]
	    pfalse[j]=f$P0[j]
	}
    }
    info$fCS[i]=log(cumprod(ptrue/pfalse))[nrow(f)]
  }
  
  if(file.exists(info[i,]$Mfile)){
    read.table(info[i,]$Mfile,header=T,stringsAsFactors=F)->m
    info$mCount[i]=nrow(m)
    ptrue=m$P0
    pfalse=m$P1
    for(j in 1:nrow(m)){
	if((!is.null(m$fetal_M[j])&&m$fetal_M[j]==2)||(is.null(m$fetal_M[j])&&m$P0<0.5)){
	    ptrue[j]=m$P1[j]
	    pfalse[j]=m$P0[j]
	}
    }
    info$mCS[i]=log(cumprod(ptrue/pfalse))[nrow(m)]
  }
}

write.table(info,"cs.path.tsv",row.names = F,quote = F,sep="\t")

