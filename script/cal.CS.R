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
    info$fCS[i]=log(cumprod(f$P0/f$P1))[nrow(f)]
  }
  
  if(file.exists(info[i,]$Mfile)){
    read.table(info[i,]$Mfile,header=T,stringsAsFactors=F)->m
    info$mCount[i]=nrow(m)
    info$mCS[i]=log(cumprod(m$P0/m$P1))[nrow(m)]
  }
}

write.table(info,"cs.tsv",row.names = F,quote = F,sep="\t")

