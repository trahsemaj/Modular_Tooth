library(ggplot2)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(FactoMineR)
#install.packages('FactoMineR')
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
select_norm_exp_upper <- function(nameList,exp){
  select_exp <- exp[toupper(rownames(exp)) %in% nameList,]
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

sample_table <- read.table('/home/james/Dropbox/Miller/data/RNA_seq_exp/total_cor/All_RNA_samples.tsv', header=TRUE,sep='\t',row.names=1)

RNA_seq_exp = read.table('/home/james/Dropbox/Miller/data/RNA_seq_exp/total_cor/kallisto_collapsed_all_genes.tsv',
                         header=TRUE,sep='\t',row.names=1)
RNA_seq_exp_ids = read.table('/home/james/Dropbox/Miller/data/RNA_seq_exp/total_cor/kallisto_collapsed_all_ensid.tsv',
                             header=TRUE,sep='\t',row.names=1)
apply(sample_table,1,function(x) {!(x %in% c('PxR_F1','CxR_F'))})
conditions_to_filter <- c('PxR_F1','CxR_F1',
                          'PAXB_DMSO','RABS_DMSO','PAXB_LDN','RABS_LDN')

samples_to_use <- rownames(sample_table)[apply(sample_table,1,function(x) {!(x[1] %in% conditions_to_filter)})]
# %in% samples_to_use
samples_to_use
colnames(RNA_seq_exp) %in% samples_to_use
RNA_seq_exp_filtered <-RNA_seq_exp[,samples_to_use]

filtered_seq <- select_norm_exp(as.vector(rownames(RNA_seq_exp_filtered)),RNA_seq_exp_filtered)
filtered_cor <- cor(filtered_seq,method='pearson')
diag(filtered_cor) <- 0

HFSC <- read.table('/home/james/Dropbox/Miller/data/BED/HFSCGS_snames.l',sep='\t')
HFSC_names = as.vector(HFSC[,c(1)])
biteCode <- read.table('/home/james/Dropbox/Miller/data/BED/BITECODE.l',sep='\t')
biteCode_names = as.vector(biteCode[,c(1)])
BC_exp <- select_norm_exp(biteCode_names,RNA_seq_exp_filtered)
BC_cor <- cor(BC_exp,method='pearson')
diag(BC_cor) <- 0
tetog <- read.table('/home/james/Dropbox/Miller/data/BED/TeToG_snames.l',sep='\t')
tetog_names = as.vector(tetog[,c(1)])
tetog_exp <- select_norm_exp_upper(c(tetog_names,'BMP6','PITX2 (2 OF 2)'),RNA_seq_exp_filtered)
lc_tetog_names <- colnames(tetog_exp)

secreted_ligands = read.table('/home/james/Dropbox/Miller/data/Networks/Ligands_v1.tsv',header=TRUE,sep='\t')
lig_names = as.vector(secreted_ligands[,c(2)])
receptors =  read.table('/home/james/Dropbox/Miller/data/Networks/Receptors_v1.tsv',header=TRUE,sep='\t')
receptors_names = as.vector(receptors[,c(2)])
TFs =  read.table('/home/james/Dropbox/Miller/data/Networks/TFactors_v1.tsv',header=TRUE,sep='\t')
TF_names = as.vector(TFs[,c(2)])
negReg =  read.table('/home/james/Dropbox/Miller/data/Networks/SecretedInhibitors_v1.tsv',header=TRUE,sep='\t')
negReg_names = as.vector(negReg[,c(2)])




###Let's find the most toothy looking genes here
gset_exp <- select_norm_exp(c(biteCode_names),RNA_seq_exp_filtered)
name_id <- data.frame(rownames(RNA_seq_exp_ids))
rownames(name_id) <- rownames(RNA_seq_exp)
name_id['bmp6',]
##changes ensid from factor to string
name_id[] <- lapply(name_id, as.character)
name_id['bmp6',]
genes_to_disp = c('bmp6','PITX2 (2 of 2)','ZNF106 (2 of 2)','plod2','snai1a')
names_to_test <- colnames(filtered_seq)
head(names_to_test)
BC_pvals <- c()
BC_diffs <- c()
for (name in names_to_test){
  filtered_cor_gset <- filtered_cor[c(colnames(gset_exp),name),c(colnames(gset_exp),name)]
  diag(filtered_cor_gset) <- 0
  pval <- wilcox.test(filtered_cor[,name],filtered_cor_gset[,name])$p.value 
  pval = log10(pval)
  BC_pvals <- c(BC_pvals,pval)
  diff <- median(filtered_cor_gset[,name])-median(filtered_cor[,name])
  BC_diffs <- c(BC_diffs,diff)
}
filtered_cor_gset <- filtered_cor[c(colnames(gset_exp),'sox2','LGR6','rspo1','rspo3','axin2','nog2'),c(colnames(gset_exp),'sox2','LGR6','rspo1','rspo3','axin2','nog2')]
diag(filtered_cor_gset) <- 0
vioplot(
  filtered_cor[,'bmp6'],filtered_cor_gset[,'bmp6'],
  filtered_cor[,'LGR6'],filtered_cor_gset[,'LGR6'],
  filtered_cor[,'rspo1'],filtered_cor_gset[,'rspo1'],
  filtered_cor[,'rspo3'],filtered_cor_gset[,'rspo3'],
  filtered_cor[,'axin2'],filtered_cor_gset[,'axin2'],
  filtered_cor[,'nog2'],filtered_cor_gset[,'nog2'],
  filtered_cor[,'sox2'],filtered_cor_gset[,'sox2'],
  drawRect = TRUE,names=c('Bmp6 All','Bmp6 BC',
                          'Lgr6 All','Lgr6 BC',
                          'Rspo1 All','Rspo1 BC',
                          'Rspo3 All','Rspo3 BC',
                          'axin2 All','axin2 BC',
                          'nog2 All','nog2 BC',
                          'sox2 All','sox2 BC'),col=c('green'))
wilcox.test( filtered_cor[,'bmp6'],filtered_cor_gset[,'bmp6'])$p.value
wilcox.test( filtered_cor[,'LGR6'],filtered_cor_gset[,'LGR6'])$p.value
wilcox.test( filtered_cor[,'rspo1'],filtered_cor_gset[,'rspo1'])$p.value
wilcox.test( filtered_cor[,'rspo3'],filtered_cor_gset[,'rspo3'])$p.value
wilcox.test( filtered_cor[,'axin2'],filtered_cor_gset[,'axin2'])$p.value
wilcox.test( filtered_cor[,'nog2'],filtered_cor_gset[,'nog2'])$p.value
wilcox.test( filtered_cor[,'sox2'],filtered_cor_gset[,'sox2'])$p.value
gset_exp <- select_norm_exp(c(HFSC_names),RNA_seq_exp_filtered)
HFSC_pvals <- c()
HFSC_diffs <- c()
for (name in names_to_test){
  filtered_cor_gset <- filtered_cor[c(colnames(gset_exp),name),c(colnames(gset_exp),name)]
  diag(filtered_cor_gset) <- 0
  pval <- wilcox.test(filtered_cor[,name],filtered_cor_gset[,name])$p.value 
  pval = log10(pval)
  HFSC_pvals <- c(HFSC_pvals,pval)
  diff <- median(filtered_cor_gset[,name])-median(filtered_cor[,name])
  HFSC_diffs <- c(HFSC_diffs,diff)
}
HFSC_diffs
HFSC_pvals
BC_pvals 
BC_diffs
pvals <- HFSC_pvals
diffs <- HFSC_diffs
sorted_names <- names_to_test[order(pvals)]
sorted_diffs <- diffs[order(pvals)]
sorted_pvals <- pvals[order(pvals)]
filtered_mean_exp = rowMeans(RNA_seq_exp_filtered)
sorted_mean_exp <- c()
sorted_ensids <- c()
for (name in sorted_names){
  cur_mean_exp <- filtered_mean_exp[name]
  sorted_mean_exp <- c(sorted_mean_exp,cur_mean_exp)
  cur_ensid <- name_id[name,]
  sorted_ensids <- c(sorted_ensids,cur_ensid)
}
length(sorted_pvals)
tail(head(sorted_pvals,n=1432))
##it looks like pval (very weakly) anticorrelates with expression
cor_df <- data.frame(sorted_ensids,sorted_names,sorted_pvals,sorted_diffs,sorted_mean_exp)
vioplot(pvals,drawRect = TRUE,col='red')
vioplot(diffs,drawRect = TRUE,col='red')

bc_diffs <- sorted_diffs[sorted_names %in% biteCode_names]
bc_pvals <- sorted_pvals[sorted_names %in% biteCode_names]
vioplot(pvals,bc_pvals,drawRect = TRUE,col='red')
vioplot(diffs,bc_diffs,drawRect = TRUE,col='red')
wilcox.test(pvals,bc_pvals)
HFSC_names[HFSC_names %in% sorted_names]
biteCode_names[biteCode_names %in% sorted_names]
sc_diffs <- sorted_diffs[sorted_names %in% HFSC_names]
sc_pvals <- sorted_pvals[sorted_names %in% HFSC_names]
vioplot(pvals,sc_pvals,drawRect = TRUE,col='red')
vioplot(diffs,sc_diffs,drawRect = TRUE,col='red')
wilcox.test(pvals,sc_pvals)


tetog_diffs <- sorted_diffs[sorted_names %in% lc_tetog_names]
tetog_pvals <- sorted_pvals[sorted_names %in% lc_tetog_names]
vioplot(pvals,tetog_pvals,drawRect = TRUE,col='red')
vioplot(diffs,tetog_diffs,drawRect = TRUE,col='red')
wilcox.test(pvals,tetog_pvals)

lig_names[lig_names %in% sorted_names]
lig_diffs <- sorted_diffs[sorted_names %in% lig_names]
lig_pvals <- sorted_pvals[sorted_names %in% lig_names]
vioplot(pvals,lig_pvals,drawRect = TRUE,col='red')
vioplot(diffs,lig_diffs,drawRect = TRUE,col='red')
wilcox.test(pvals,lig_pvals)

receptors_names[receptors_names %in% sorted_names]
receptors_diffs <- sorted_diffs[sorted_names %in% receptors_names]
receptors_pvals <- sorted_pvals[sorted_names %in% receptors_names]
vioplot(pvals,receptors_pvals,drawRect = TRUE,col='red')
vioplot(diffs,receptors_diffs,drawRect = TRUE,col='red')
wilcox.test(pvals,receptors_pvals)

TF_names[TF_names %in% sorted_names]
TF_diffs <- sorted_diffs[sorted_names %in% TF_names]
TF_pvals <- sorted_pvals[sorted_names %in% TF_names]
vioplot(pvals,TF_pvals,drawRect = TRUE,col='red')
vioplot(diffs,TF_diffs,drawRect = TRUE,col='red')
wilcox.test(pvals,TF_pvals)

negReg_diffs <- sorted_diffs[sorted_names %in% negReg_names]
negReg_pvals <- sorted_pvals[sorted_names %in% negReg_names]
vioplot(pvals,negReg_pvals,drawRect = TRUE,col='red')
vioplot(diffs,negReg_diffs,drawRect = TRUE,col='red')
wilcox.test(pvals,negReg_pvals)

wilcox.test(lig_pvals,receptors_pvals)
wilcox.test(lig_pvals,TF_pvals)



vioplot(pvals,bc_pvals,tetog_pvals,sc_pvals,drawRect = TRUE,
        names = c('all','BiteCode','TeTog','HFSC'),col='blue')
wilcox.test(pvals,sc_pvals)
vioplot(pvals,tetog_pvals,drawRect = TRUE,
        names = c('all','TeTog'),col='red')

vioplot(pvals,TF_pvals,receptors_pvals,lig_pvals,negReg_pvals,drawRect = TRUE,
        names = c('all','TFs','Receptors','Ligand','NegReg'),col='blue')


negReg_pvals
vioplot(diffs,bc_diffs,drawRect = TRUE,col='red')
plot(log(sorted_mean_exp)~sorted_pvals,data=cor_df)
summary(lm(log(sorted_mean_exp)~sorted_pvals,data=cor_df))
abline(lm(log(sorted_mean_exp)~sorted_pvals,data=cor_df))
filtered_mean_exp[1]


filtered_cor_gset <- filtered_cor[c(colnames(gset_exp),'sox2','dcn'),c(colnames(gset_exp),'sox2','dcn')]
vioplot(
  filtered_cor[,'sox2'],filtered_cor_gset[,'sox2'],
  filtered_cor[,'bmp6'],filtered_cor_gset[,'bmp6'],
  filtered_cor[,'dcn'],filtered_cor_gset[,'dcn'],
  drawRect = TRUE,names=c('Sox2 All','Sox2 HFSC',
                           'Bmp6 All','Bmp6 HFSC',
                          'Dcn All','Dcn HFSC'),col=c('green'))
###make specific plots hers
write.table(sorted_names,file='/home/james/Dropbox/Miller/data/RNA_seq_exp/total_cor/HFSC_names_allsamples.l',row.names = FALSE,col.names = FALSE,quote=FALSE)
write.table(cor_df,file='/home/james/Dropbox/Miller/data/RNA_seq_exp/total_cor/BC_filtered_sorted_genes_pvals_diffs_032017.tsv',row.names = FALSE,col.names = TRUE,quote=FALSE,sep='\t')

