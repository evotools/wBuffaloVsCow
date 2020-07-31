library(ggplot2)
library(ggpubr)
library(RColorBrewer)

dt=read.table("metadata.txt",sep = "\t",header = T)
PCA=read.table("PCA_new_20.eigenvec",sep = "\t",header = T)
new_dt=merge(dt,PCA,by="FID")
new_dt$IID=NULL


#################################For creating png plot############################
customPalette=c("#E41A1C","#377EB8","#4DAF4A","#A65628","#984EA3","#FF7F00","#FFFF33")
png('PCA_new_20_with_med_bhad_pand_bigger_fonts.png', units="in", width=18, height=10, res=300)
ggplot(new_dt, aes(PC1, PC2,color=Breed,shape=Sequencing.centre)) +theme_pubclean()+
  geom_point(size=5) + labs(x="PC1",y="PC2")+
  theme(text = element_text(size=22))+labs(shape="Sequencing centre")+
  scale_color_manual(values=customPalette)
dev.off()
