library(FactoMineR)
#install.packages('ggrepel')
library(ggplot2)
library(ggrepel)

Dropbox = 'C:/Users/trahs/Dropbox/'

MOAV_gts = read.table(paste(Dropbox,'Miller/data/Genomic_Analysis/MOAV_071817_tested_noFRI_filtered_haplotype.gt',sep=""),header=TRUE,sep='\t')
MOAV_gts <- MOAV_gts[,c(-1)]

#I4_gts = read.table('/home/james/Data/MOAV_vcfs/MOAV_071817_tested_filtered_haplotype_I4.gt',header=TRUE,sep='\t')
#I4_gts <- I4_gts[,c(-1)]

M <- data.matrix(MOAV_gts)
pca_res <- PCA(t(M),scale.unit=FALSE,ncp=10,graph=FALSE)

pca_df <- data.frame(pca_res$ind$coord)

colnames(M)

tested_pop_names = c("CERC.wd454", "BEPA.wdBxL", "LITC.wmBxL", "FTC.wdFTxL", "LITC.mFTxL", "LITC.wmL28",
                    "LITC.wmL29", "PAXB.wdPxL", "FTC.lmBRev", "LITC.ldBRv", "PAXB.ldPxR", "PAXB.dPxJa", "PAXB.dPxJl",
                    "PAXB.dWVir", "CERC.dCxRB", "JAMA.mPxJl", "Hutu.F3",    "CCD.2007",   "LITC.wmPB",  "EnosB.wdPB",
                    "PriestB.wd", "PaxtonB.wd", "RABS.lmCxR")

tested_pop_names = c("CERC1", "BEPA", "LITC6", "FTC1", "LITC4", "LITC1",
                     "LITC2", "PAXB1", "FTC2", "LITC5", "PAXB2", "PAXB3", "PAXB4",
                     "PAXB5", "CERC2", "JAMA", "HUTU",    "CCD",   "LITC3",  "ENOB",
                     "PRIB", "PAXB6", "RABS")


tested_pop_colors= c("darkgreen", "darkblue", "red", "deepskyblue3", "red", "red",
                     "red", "blue", "deepskyblue3", "red", "blue", "blue", "blue",
                     "blue", "darkgreen", "darkred", "darkcyan",    "darkcyan",   "red",  "lightslateblue",
                     "lightslateblue", "lightslateblue", "firebrick1")

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC2 ",round(pca_res$eig[,"percentage of variance"][2],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.2)) +
  geom_point(color = tested_pop_colors) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
ggsave(paste(Dropbox,'Miller/figures/Multipop_intron4/Genomic_PCA_1_2.pdf',sep=""),width=6,height=4.5,units='in')


xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC3 ",round(pca_res$eig[,"percentage of variance"][3],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.3)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_point(color = tested_pop_colors) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
ggsave(paste(Dropbox,'Miller/figures/Multipop_intron4/Genomic_PCA_1_3.pdf',sep=""),width=6,height=4.5,units='in')


xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC4 ",round(pca_res$eig[,"percentage of variance"][4],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.4)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_point(color = tested_pop_colors) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
ggsave(paste(Dropbox,'Miller/figures/Multipop_intron4/Genomic_PCA_1_4.pdf',sep=""),width=6,height=4.5,units='in')

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC5 ",round(pca_res$eig[,"percentage of variance"][5],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.5)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_point(color = tested_pop_colors) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
ggsave(paste(Dropbox,'Miller/figures/Multipop_intron4/Genomic_PCA_1_5.pdf',sep=""),width=6,height=4.5,units='in')

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC6 ",round(pca_res$eig[,"percentage of variance"][6],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.6)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_point(color = tested_pop_colors) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16)

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC7 ",round(pca_res$eig[,"percentage of variance"][7],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.7)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_point(color = tested_pop_colors) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16)

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC8 ",round(pca_res$eig[,"percentage of variance"][8],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.8)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_point(color = tested_pop_colors) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16)

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC9 ",round(pca_res$eig[,"percentage of variance"][9],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.9)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_point(color = tested_pop_colors) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16)

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC10 ",round(pca_res$eig[,"percentage of variance"][10],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.10)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_point(color = tested_pop_colors) +
  geom_text_repel(aes(label = tested_pop_names)) +
  theme_classic(base_size = 16)





xlabel = paste("PC1 ",round(pca_res$eig$"percentage of variance"[1],digits=2),"%")
ylabel = paste("PC2 ",round(pca_res$eig$"percentage of variance"[2],digits=2),"%")
plot(pca_res$ind$coord[,c('Dim.1')],pca_res$ind$coord[,c('Dim.2')],ylab=ylabel,xlab=xlabel,pch=19,xlim=c(-250,-150),ylim=c(-100,-50))
text(pca_res$ind$coord[,c('Dim.1')],pca_res$ind$coord[,c('Dim.2')],labels=colnames(M),pos=3,offset=.4,cex=.7)

xlabel = paste("PC1 ",round(pca_res$eig$"percentage of variance"[1],digits=2),"%")
ylabel = paste("PC3 ",round(pca_res$eig$"percentage of variance"[3],digits=2),"%")
plot(pca_res$ind$coord[,c('Dim.1')],pca_res$ind$coord[,c('Dim.3')],ylab=ylabel,xlab=xlabel,pch=19)
text(pca_res$ind$coord[,c('Dim.1')],pca_res$ind$coord[,c('Dim.3')],labels=colnames(M),pos=3,offset=.4,cex=.7)
