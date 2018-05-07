library(ggplot2)
library(ggplot2)
library(RColorBrewer)
#install.packages('vioplot')
library(FactoMineR)
library(ggrepel)
library(vioplot)
my_palette <- colorRampPalette(c("green",'black','red'))(n = 100)
select_norm_exp <- function(nameList,exp){
  select_exp <- exp[rownames(exp) %in% nameList,]
  #select_exp <- select_exp[apply(select_exp, 1, function(x) length(x[x>.3])>=(length(colnames(exp)))-5),]
  select_exp  <- t(as.matrix(select_exp ))
  #select_exp  <- log2(select_exp )
  ##select_exp  = scale(select_exp ,scale=TRUE,center=TRUE)
  ##select_exp  = scale(select_exp ,scale=TRUE,center=FALSE)
  select_exp  = scale(select_exp ,scale=colSums(select_exp),center=FALSE)
  select_exp <- select_exp[,colSums(is.na(select_exp)) < 1]
  ##select_exp  = scale(select_exp ,scale=colSums(select_exp),center=TRUE)
  return(select_exp)
}
Dropbox = 'C:/Users/trahs/Dropbox/'
paste(Dropbox,'',sep="")
sample_table <- read.table(paste(Dropbox,'Miller/data/RNA_seq_exp/total_cor/All_RNA_samples.tsv',sep=""), header=TRUE,sep='\t',row.names=1)

RNA_seq_exp = read.table(paste(Dropbox,'Miller/data/RNA_seq_exp/total_cor/kallisto_collapsed_all_genes.tsv',sep=""),
                         header=TRUE,sep='\t',row.names=1)
RNA_seq_exp_ids = read.table(paste(Dropbox,'Miller/data/RNA_seq_exp/total_cor/kallisto_collapsed_all_ensid.tsv',sep=""),
                             header=TRUE,sep='\t',row.names=1)
apply(sample_table,1,function(x) {!(x %in% c('PxR_F1','CxR_F'))})
conditions_to_filter <- c('PxR_F1','CxR_F1',
                          'PAXB_DMSO','RABS_DMSO','PAXB_LDN','RABS_LDN')

samples_to_use <- rownames(sample_table)[apply(sample_table,1,function(x) {!(x[1] %in% conditions_to_filter)})]
# %in% samples_to_use
samples_to_use
colnames(RNA_seq_exp) %in% samples_to_use
RNA_seq_exp_filtered <-RNA_seq_exp[,samples_to_use]



RNA_samples <- colnames(RNA_seq_exp)
RNA_samples_filtered <- colnames(RNA_seq_exp_filtered)
RNA_samples_filtered 
RNA_samples_colors <- c("darkseagreen3","darkseagreen3","darkseagreen3",
                        "red","red","red",
                        "blue","blue","blue",
                        "darkblue","darkblue","darkblue",
                        "darkred","darkred","darkred",
                        "midnightblue","midnightblue","midnightblue",
                        "orangered","orangered","orangered",
                        "lightseagreen","lightseagreen","lightseagreen",
                        "midnightblue","midnightblue","midnightblue",
                        "blue","blue","blue",
                        "red4","red4","red4",
                        "midnightblue","midnightblue","midnightblue",
                        "midnightblue","midnightblue","midnightblue",
                        "blue","blue","blue",
                        "purple","purple","purple","purple","purple",
                        "orange","orange","orange","orange","orange",
                        "lightcoral","lightcoral","lightcoral"
)
RNA_samples_colors_filtered <- c("darkseagreen3","darkseagreen3","darkseagreen3",
                                 "red","red","red",
                                 "blue","blue","blue",
                                 "lightseagreen","lightseagreen","lightseagreen",
                                 "blue","blue","blue",
                                 "blue","blue","blue",
                                 "red4","red4","red4",
                                 "blue","blue","blue",
                                 "blue","blue","blue",
                                 "blue","blue","blue",
                                 "lightcoral","lightcoral","lightcoral"
)
RNA_samples_shape_filtered <- c(16,16,16,
                                16,16,16,
                                16,16,16,
                                16,16,16,
                                2,2,2,
                                17,17,17,
                                16,16,16,
                                0,0,0,
                                7,7,7,
                                15,15,15,
                                16,16,16

                                
)
colnames(M)
M = data.matrix(RNA_seq_exp)
M = M + .1
M  = scale(M ,scale=TRUE,center=FALSE)
M <-log2(M)
pca_res <- PCA(t(M),scale.unit=FALSE,ncp=10,graph=FALSE)
pca_df <- data.frame(pca_res$ind$coord)
pca_res$eig[,"percentage of variance"]


xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC2 ",round(pca_res$eig[,"percentage of variance"][2],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.2)) +
  geom_point(color = RNA_samples_colors) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/RNA_seq_PCA_1_2.pdf',width=6,height=4.5,units='in')


xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC3 ",round(pca_res$eig[,"percentage of variance"][3],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.3)) +
  geom_point(color = RNA_samples_colors) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/RNA_seq_PCA_1_3.pdf',width=6,height=4.5,units='in')


xlabel = paste("PC2 ",round(pca_res$eig[,"percentage of variance"][2],digits=2),"%")
ylabel = paste("PC3 ",round(pca_res$eig[,"percentage of variance"][3],digits=2),"%")
ggplot(pca_df,aes(Dim.2,Dim.3)) +
  geom_point(color = RNA_samples_colors) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/RNA_seq_PCA_2_3.pdf',width=6,height=4.5,units='in')

M = data.matrix(RNA_seq_exp_filtered)
M = M + .1
M  = scale(M ,scale=TRUE,center=FALSE)
M <-log2(M)
pca_res <- PCA(t(M),scale.unit=FALSE,ncp=10,graph=FALSE)
pca_df <- data.frame(pca_res$ind$coord)
pca_res$eig[,"percentage of variance"]


xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC2 ",round(pca_res$eig[,"percentage of variance"][2],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.2)) +
  geom_point(color = RNA_samples_colors_filtered,shape=RNA_samples_shape_filtered,size=4) +
  xlab(xlabel) + 
  ylab(ylabel) +
  #geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)
#ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/RNA_seq_filtered_PCA_1_2.pdf',width=6,height=4.5,units='in')
ggsave(paste(Dropbox,'Miller/figures/Multipop_intron4/RNA_seq_filtered_PCA_1_2_noNames.pdf',sep=""),width=6,height=4.5,units='in')


xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC3 ",round(pca_res$eig[,"percentage of variance"][3],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.3)) +
  geom_point(color = RNA_samples_colors_filtered,shape=RNA_samples_shape_filtered,size=4) +
  xlab(xlabel) + 
  ylab(ylabel) +
  #geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)
#ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/RNA_seq_filtered_PCA_1_3.pdf',width=6,height=4.5,units='in')
ggsave(paste(Dropbox,'Miller/figures/Multipop_intron4/RNA_seq_filtered_PCA_1_3_noNames.pdf',sep=""),width=6,height=4.5,units='in')

xlabel = paste("PC2 ",round(pca_res$eig[,"percentage of variance"][2],digits=2),"%")
ylabel = paste("PC3 ",round(pca_res$eig[,"percentage of variance"][3],digits=2),"%")
ggplot(pca_df,aes(Dim.2,Dim.3)) +
  geom_point(color = RNA_samples_colors_filtered) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/RNA_seq_filtered_PCA_2_3.pdf',width=6,height=4.5,units='in')

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC4 ",round(pca_res$eig[,"percentage of variance"][4],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.4)) +
  geom_point(color = RNA_samples_colors_filtered) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)

xlabel = paste("PC2 ",round(pca_res$eig[,"percentage of variance"][2],digits=2),"%")
ylabel = paste("PC4 ",round(pca_res$eig[,"percentage of variance"][4],digits=2),"%")
ggplot(pca_df,aes(Dim.2,Dim.4)) +
  geom_point(color = RNA_samples_colors_filtered) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)

xlabel = paste("PC3 ",round(pca_res$eig[,"percentage of variance"][3],digits=2),"%")
ylabel = paste("PC4 ",round(pca_res$eig[,"percentage of variance"][4],digits=2),"%")
ggplot(pca_df,aes(Dim.3,Dim.4)) +
  geom_point(color = RNA_samples_colors_filtered) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)

xlabel = paste("PC1 ",round(pca_res$eig[,"percentage of variance"][1],digits=2),"%")
ylabel = paste("PC5 ",round(pca_res$eig[,"percentage of variance"][5],digits=2),"%")
ggplot(pca_df,aes(Dim.1,Dim.5)) +
  geom_point(color = RNA_samples_colors_filtered) +
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_text_repel(aes(label = colnames(M))) +
  theme_bw(base_size = 16)