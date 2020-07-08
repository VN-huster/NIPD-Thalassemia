scan("snp.count.txt")->snp
data.frame(snp)->snp
#plot(density(snp$snp),col=3,cex=2,lwd=3,xlab="SNP count",ylab="Density");abline(v=20)
library(ggplot2)
ggplot(snp,aes(snp))+
  xlab("SNP count")+
  ylab("Density")+
  geom_density(col=3,size=2)+
  geom_vline(xintercept = 20)+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.title = element_text(face = "bold"),
    text = element_text(size=20)
  )
ggsave("snp.count.tiff",dpi=600,device = "tiff",compression = "lzw")
ggsave("snp.count.png",dpi=600,device = "png",compression = "lzw")
