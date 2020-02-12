# HBA
gene.chr="chr16"
gene.start=c(222846)
gene.end=c(227520)
gene.name="HBA1-HBA2"

extL=2e5
xl=c(0,4e5)

col.unmatch="grey"
col.matchF="lightseagreen"
col.matchM="lightcoral"

pdf("FigureS2_alpha.v6.pdf",width=15,height=25)
#tiff("alpha_haplotying_consistency.v2.tiff",width=10,height=14,units='in',res=600,compression='lzw')
#jpeg(filename ="alpha_haplotying_consistency_arial.v1_nomargin.jpeg",width=10,height=14,units='in',res=600,family='Arial',pointsize=8)
par(mar=c(5,2,0,0),mgp=c(3,.5,0),mfrow=c(1,2))
read.table("alpha_file.list",stringsAsFactors = F)->flist

#mat
plot(0,t='n',
     xlab="Position (chromosome 16)/Mbp",
     ylab="",
     cex.lab=1.5,
     xlim=xl,ylim=c(1,nrow(flist)+4.5),
     bty='n',
     xaxt='n',yaxt='n',
     xpd=T)
xls=seq(0,xl[2],by=1e5)
axis(1,at=xls,labels = xls/1e6,cex.axis=1.3)
for(i in 1:nrow(flist)){
    if ( i %% 2 == 0){
	read.table(flist[i,2],stringsAsFactors = T,header=T)->t
	t$match=col.unmatch
	nrow(t)
	t[t$parent==t$fetus,]$match=col.matchM
	t[t$pos<=gene.end[1]+extL&t$pos>=0,]->t
	t$start=0
	t$end=gene.end[1]+extL
	points(t$pos,rep(nrow(flist)+1-i+1-0.1,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
	polygon(x=c(xl[1]+5e4,xl[2]-5e4),y=c(nrow(flist)+1-i+1-0.1,nrow(flist)+1-i+1-0.1),col='grey65',density=-1,border='grey65',lwd=.5)
	#text(3e4,i+0.2-1,labels =paste0(flist[nrow(flist)+1-i,1]),pos=2,cex=1.5)
	points(t$pos,rep(nrow(flist)+1-i+0.3,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
	polygon(x=c(xl[1]+5e4,xl[2]-5e4),y=c(nrow(flist)+1-i+0.3,nrow(flist)+1-i+0.3),col='grey65',density=-1,border='grey65',lwd=.5)
	text(3e4,i+0.5-1,labels =paste0(flist[nrow(flist)+1-i,1]),pos=2,cex=1.5)
    }
}
rect(gene.start,0.5,gene.end,nrow(flist)+0.5,col= NA, border='lightslateblue')
text(mean(gene.start,gene.end),nrow(flist)+1,labels = gene.name,cex=1.5)
text(1.9e5,nrow(flist)+2,labels = "Maternal haplotype",cex=2,font=2)

#pat
plot(0,t='n',
     xlab="Position (chromosome 16)/Mbp",
     ylab="",
     cex.lab=1.5,
     xlim=xl,ylim=c(1,nrow(flist)+4.5),
     bty='n',
     xaxt='n',yaxt='n',
     xpd=T)
xls=seq(0,xl[2],by=1e5)
axis(1,at=xls,labels = xls/1e6,cex.axis=1.3)

for(i in 1:nrow(flist)){
    if ( i %% 2 != 0){
	read.table(flist[i,2],stringsAsFactors = T,header=T)->t
	t$match=col.unmatch
	nrow(t)
	t[t$parent==t$fetus,]$match=col.matchF
	t[t$pos<=gene.end[1]+extL&t$pos>=0,]->t
	t$start=0
	t$end=gene.end[1]+extL
	points(t$pos,rep(nrow(flist)+1-i-0.1,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
	polygon(x=c(xl[1]+5e4,xl[2]-5e4),y=c(nrow(flist)+1-i-0.1,nrow(flist)+1-i-0.1),col='grey65',density=-1,border='grey65',lwd=.5)
	#text(3e4,i+0.2,labels =paste0(flist[nrow(flist)+1-i,1]),pos=2,cex=1.5)
	points(t$pos,rep(nrow(flist)-i+0.3,nrow(t)),t='p',pch=20,col=t$match,cex=1.5)
	polygon(x=c(xl[1]+5e4,xl[2]-5e4),y=c(nrow(flist)-i+0.3,nrow(flist)-i+0.3),col='grey65',density=-1,border='grey65',lwd=.5)
	text(3e4,i+0.5,labels =paste0(flist[nrow(flist)-i,1]),pos=2,cex=1.5)						    }
}
rect(gene.start,0.5,gene.end,nrow(flist)+0.5,col= NA, border='lightslateblue')
text(mean(gene.start,gene.end),nrow(flist)+1,labels = gene.name,cex=1.5)
text(1.9e5,nrow(flist)+2,labels = "Paternal haplotype",cex=2,font=2)


legend(1.8e5,nrow(flist)+6.3,legend=c("Paternal matched SNPs","Maternal matched SNPs","Parental mismatched SNPs"),pch=20,cex=1.5,col=c(col.matchF,col.matchM,col.unmatch))
dev.off()
