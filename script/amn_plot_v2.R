args=commandArgs(T)
if(length(args)<2){
    print("--args file gene [pdf:file.pdf]")
    q()
}
file=args[1]
gene=args[2]
if(length(args)>2){
    pdf=args[3]
}else{
    pdf=paste0(file,".pdf")
}
read.table(file,sep="\t",stringsAsFactors=F)->Data
#colnames(Data)=c("Chr","Pos","AMN_F","P0","P1","fetal_F","Note","Gene");
colnames(Data)=c("Chr","Pos","P0","P1","fetal_F", "amn", "Note","Gene");
#read.table("/vol6/home/bgi_wangyy/project/PGD/NIPTv2/NIPT_V2_gene.txt",sep="\t",stringsAsFactors=F)->Gene
read.table(gene,sep="\t",stringsAsFactors=F)->Gene

pdf(pdf,width=22)

for(i in 1:nrow(Gene)){
    Data[Data$Gene==Gene[i,6],]->tmp;
#    unique(tmp$Note)->Note
    paste0(Gene[i,6],c("_US_500K-1M","_US_500K","","_DS_500K","_DS_500K-1M"))->Note
    if(nrow(tmp)>0){
        par(mar=c(3,0,3,0),mfrow=c(1,length(Note)))
        for(j in 1:length(Note)){
            tmp[tmp$Note==Note[j],]->tmp2
            if(nrow(tmp2)>0){
                plot(tmp2$P0-tmp2$P1~tmp2$Pos,xlab="pos",ylab="P0-P1",main=Note[j],pch=20,cex=1,col=2,ylim=c(-1,1))
                points(1.5-tmp2$fetal_F~tmp2$Pos,col=3,t="b",cex=1,lwd=4)
                points(1.3-tmp2[tmp2$amn!=1.5,]$amn~tmp2[tmp2$amn!=1.5,]$Pos,col=4,t="b",cex=1,lwd=4)
#                points(1.3-tmp2$AMN_F~tmp2$Pos,col=4,t="b")
                abline(h=0,col=5)
            }else{
                plot(0,xlab="pos",ylab="P0-P1",main=Note[j],pch=20,cex=0.5,col=2,ylim=c(-1,1),xaxt='n',t='n')
                abline(h=0,col=5)
            }
        }
    }
}

dev.off()

# print "$Chr\t$Pos\t$AMN_F\t$AMN_M\t$P0\t$P1\t$fetal_F\t$fetal_M\t".$hash{"$ln[0]\t$ln[1]"}."\n";
