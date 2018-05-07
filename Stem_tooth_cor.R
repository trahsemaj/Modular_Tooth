library(ggplot2)
library(RColorBrewer)
library(gplots)
library(FactoMineR)
library(vioplot)
install.packages('beanplot')
library(beanplot)
source_gist()
require(digest)
require(devtools)
install.packages('devtools')
#source_gist("https://gist.github.com/mbjoseph/5852613")

my_palette <- colorRampPalette(c("green",'black','red'))(n = 100)
select_norm_exp <- function(nameList,exp){
  select_exp <- exp[rownames(exp) %in% nameList,]
  #select_exp <- select_exp[apply(select_exp, 1, function(x) length(x[x>.3])>=(length(colnames(exp)))-5),]
  select_exp  <- t(as.matrix(select_exp ))
  #select_exp  <- log2(select_exp )
  select_exp  = scale(select_exp ,scale=TRUE,center=TRUE)
  ##select_exp  = scale(select_exp ,scale=TRUE,center=FALSE)
  #select_exp  = scale(select_exp ,scale=colSums(select_exp),center=FALSE)
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
Dropbox = 'C:/Users/trahs/Dropbox/'
Dropbox = '/home/james/Dropbox/'
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

filtered_seq <- select_norm_exp(as.vector(rownames(RNA_seq_exp_filtered)),RNA_seq_exp_filtered)
filtered_cor <- cor(filtered_seq,method='pearson')
diag(filtered_cor) <- 0
head(sort(filtered_cor['PITX2 (2 of 2)',],decreasing = TRUE))

total_seq <- select_norm_exp(as.vector(rownames(RNA_seq_exp)),RNA_seq_exp)
total_cor <- cor(total_seq,method='pearson')
diag(total_cor) <- 0


HFSC <- read.table(paste(Dropbox,'Miller/data/BED/HFSCGS_snames.l',sep=""),sep='\t')
HFSC_names = as.vector(HFSC[,c(1)])
biteCode <- read.table(paste(Dropbox,'Miller/data/BED/BITECODE.l',sep=""),sep='\t')
biteCode_names = as.vector(biteCode[,c(1)])
BC_exp <- select_norm_exp(biteCode_names,RNA_seq_exp_filtered)
BC_cor <- cor(BC_exp,method='pearson')
diag(BC_cor) <- 0
tetog <- read.table(paste(Dropbox,'Miller/data/BED/TeToG_snames.l',sep=""),sep='\t')
tetog_names = as.vector(tetog[,c(1)])
tetog_exp <- select_norm_exp_upper(c(tetog_names,'BMP6','PITX2 (2 OF 2)'),RNA_seq_exp_filtered)
lc_tetog_names <- colnames(tetog_exp)

secreted_ligands = read.table(paste(Dropbox,'Miller/data/Networks/Ligands_v1.tsv',sep=""),header=TRUE,sep='\t')
lig_names = as.vector(secreted_ligands[,c(2)])
receptors =  read.table(paste(Dropbox,'Miller/data/Networks/Receptors_v1.tsv',sep=""),header=TRUE,sep='\t')
receptors_names = as.vector(receptors[,c(2)])
TFs =  read.table(paste(Dropbox,'Miller/data/Networks/TFactors_v1.tsv',sep=""),header=TRUE,sep='\t')
TF_names = as.vector(TFs[,c(2)])
negReg =  read.table(paste(Dropbox,'Miller/data/Networks/SecretedInhibitors_v1.tsv',sep=""),header=TRUE,sep='\t')
negReg_names = as.vector(negReg[,c(2)])

bone_sc =  read.table(paste(Dropbox,'Miller/data/Networks/Dental_Bone_mes_sc.tsv',sep=""),header=TRUE,sep='\t')
bone_sc_names = as.vector(bone_sc[,c(2)])

tooth_sc =  read.table(paste(Dropbox,'Miller/data/Networks/Dental_mes_sc.tsv',sep=""),header=TRUE,sep='\t')
tooth_sc_names = as.vector(tooth_sc[,c(2)])



##get BC pvals
gset_exp <- select_norm_exp(c(biteCode_names),RNA_seq_exp_filtered)
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

##get all genes (no filter)
gset_exp_all <- select_norm_exp(c(biteCode_names),RNA_seq_exp)
names_to_test <- colnames(total_seq)
head(names_to_test)
total_BC_pvals <- c()
total_BC_diffs <- c()
for (name in names_to_test){
  total_cor_gset <- total_cor[c(colnames(gset_exp_all),name),c(colnames(gset_exp_all),name)]
  diag(total_cor_gset) <- 0
  pval <- wilcox.test(total_cor[,name],total_cor_gset[,name])$p.value 
  pval = log10(pval)
  total_BC_pvals <- c(total_BC_pvals,pval)
  diff <- median(total_cor_gset[,name])-median(total_cor[,name])
  total_BC_diffs <- c(total_BC_diffs,diff)
}


##Get HFSC pvals
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


##Test HFSC correlates
gset_exp_all <- select_norm_exp(c(biteCode_names),RNA_seq_exp)
names_to_test <- colnames(total_seq)
total_pvals <- total_BC_pvals


total_sorted_names <- names_to_test[order(total_pvals)]

total_sorted_pvals <- pvals[order(total_pvals)]
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

pvals <- BC_pvals
diffs <- BC_diffs
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


lc_tetog_names
BC_unique_names <- biteCode_names[!(biteCode_names %in% HFSC_names)]
HFSC_unique_names <- HFSC_names[!(HFSC_names %in% biteCode_names)]
bc_diffs <- sorted_diffs[sorted_names %in% BC_unique_names]
bc_pvals <- sorted_pvals[sorted_names %in% BC_unique_names]
sc_diffs <- sorted_diffs[sorted_names %in% HFSC_unique_names]
sc_pvals <- sorted_pvals[sorted_names %in% HFSC_unique_names]
tetog_diffs <- sorted_diffs[sorted_names %in% lc_tetog_names]
tetog_pvals <- sorted_pvals[sorted_names %in% lc_tetog_names]
lig_diffs <- sorted_diffs[sorted_names %in% lig_names]
lig_pvals <- sorted_pvals[sorted_names %in% lig_names]
receptors_diffs <- sorted_diffs[sorted_names %in% receptors_names]
receptors_pvals <- sorted_pvals[sorted_names %in% receptors_names]
TF_diffs <- sorted_diffs[sorted_names %in% TF_names]
TF_pvals <- sorted_pvals[sorted_names %in% TF_names]
negReg_diffs <- sorted_diffs[sorted_names %in% negReg_names]
negReg_pvals <- sorted_pvals[sorted_names %in% negReg_names]
bone_sc_diffs <- sorted_diffs[sorted_names %in% bone_sc_names]
bone_sc_pvals <- sorted_pvals[sorted_names %in% bone_sc_names]
tooth_sc_diffs <- sorted_diffs[sorted_names %in% tooth_sc_names]
tooth_sc_pvals <- sorted_pvals[sorted_names %in% tooth_sc_names]

sorted_names[sorted_names %in% lc_tetog_names]
vioplot(pvals,sc_pvals,bc_pvals,tetog_pvals,drawRect = TRUE,
        names = c('all','HFSC','BiteCode','TeTog'),col='blue')

vioplot(pvals,tetog_pvals,drawRect = TRUE,
        names = c('all','TeTog'),col=c('red','blue'))

dev.off()
beanplot(pvals,tetog_pvals,side='both',
        names = c('all','TeTog'),col = list("lightblue", c("red", "black")),ll=.01)
dev.print(pdf, paste(Dropbox,'Miller/figures/Multipop_intron4/Genome_wide_TeTog_splitvio.pdf',sep=""))

gw_d <- density(pvals)
plot(gw_d, type="n", main="Filtered",xlim=c(0,-35),xlab='Log10 (Mann-Whitney U p-vals)')
polygon(gw_d, col=rgb(0,0,1,.5), border="black")
#rug(pvals, col="blue")

tetog_d <- density(tetog_pvals)
#plot(gw_d, type="n", main="test",xlim=c(0,-35))
polygon(tetog_d, col=rgb(1,0,0,.5), border="black")
rug(tetog_pvals, col="red")
dev.print(pdf, paste(Dropbox,'Miller/figures/Multipop_intron4/Genome_wide_TeTog_density.pdf',sep=""))

##total_cor

total_tetog_pvals <- total_sorted_pvals[total_sorted_names %in% lc_tetog_names]
gw_d <- density(total_pvals)
plot(gw_d, type="n", main="Total",xlim=c(0,-35),xlab='Log10 (Mann-Whitney U p-vals)')
polygon(gw_d, col=rgb(0,0,1,.5), border="black")
#rug(pvals, col="blue")

tetog_d <- density(total_tetog_pvals)
#plot(gw_d, type="n", main="test",xlim=c(0,-35))
polygon(tetog_d, col=rgb(1,0,0,.5), border="black")
rug(tetog_pvals, col="red")


vioplot(filtered_cor[,'bmp6'],range=.1,h=.05,drawRect = TRUE,ylim=c(-.8,1),names=c('Bmp6'),col=c('blue'))
title(ylab='CorCoef',ps=20)
Pitx2_cor = filtered_cor['PITX2 (2 of 2)','bmp6'] 
Plod2_cor = filtered_cor['plod2','bmp6'] 
Lgr6_cor = filtered_cor['LGR6','bmp6'] 
Bmp2b_cor = filtered_cor['bmp2b','bmp6'] 
Bmp2b_cor
Plod2_cor
Pitx2_cor
Lgr6_cor

filtered_cor['calb2a','bmp6'] 
total_cor['calb2a','bmp6'] 

arrows(.62,1,x1=.99,y1=Bmp2b_cor)
text(.55,1,labels='Bmp2b',font=3)

arrows(.62,.9,x1=.99,y1=Plod2_cor)
text(.55,.9,labels='Plod2',font=3)

arrows(.62,.8,x1=.99,y1=Pitx2_cor)
text(.55,.8,labels='Pitx2',font=3)

arrows(.62,.7,x1=.99,y1=Lgr6_cor)
text(.55,.7,labels='Lgr6',font=3)

filtered_cor['sox2','bmp6'] 
filtered_cor['spp1','bmp6'] 
tail(head(sort(filtered_cor[,'bmp6'],decreasing = TRUE),n=2500) )

vioplot(total_cor[,'bmp6'],range=.1,h=.05,drawRect = TRUE,ylim=c(-.8,1),names=c('Bmp6'),col=c('blue'))
title(ylab='CorCoef',ps=20)
Pitx2_cor = total_cor['PITX2 (2 of 2)','bmp6'] 
Plod2_cor = total_cor['plod2','bmp6'] 
Lgr6_cor = total_cor['LGR6','bmp6'] 
Bmp4_cor = total_cor['bmp4','bmp6'] 
Bmp2b_cor = total_cor['bmp2b','bmp6'] 

Bmp4_cor
Plod2_cor
Pitx2_cor
Lgr6_cor
total_cor['sox2','bmp6'] 
total_cor['spp1','bmp6'] 
arrows(.62,1,x1=.99,y1=Lgr6_cor)
text(.55,1,labels='Lgr6',font=3)

arrows(.62,.9,x1=.99,y1=Bmp2b_cor)
text(.55,.9,labels='Bmp2b',font=3)

arrows(.62,.8,x1=.99,y1=Plod2_cor)
text(.55,.8,labels='Plod2',font=3)

arrows(.62,.7,x1=.99,y1=Pitx2_cor)
text(.55,.7,labels='Pitx2',font=3)

head(sort(total_cor[,'bmp6'],decreasing = TRUE),n=50) 


dev.print(pdf, paste(Dropbox,'Miller/figures/Multipop_intron4/HFSC_BC_tetog_pvals_vio.pdf',sep=""))
vioplot(pvals,TF_pvals,receptors_pvals,lig_pvals,negReg_pvals,drawRect = TRUE,
        names = c('all','TFs','Receptors','Ligand','NegReg'),col='blue')
dev.print(pdf, paste(Dropbox,'Miller/figures/Multipop_intron4/HFSC_TF_ligands_pvals_vio.pdf',sep=""))


wilcox.test(pvals,bc_pvals)        ##0.0001372
wilcox.test(pvals,tetog_pvals)     ##0.3678
wilcox.test(pvals,sc_pvals)        ##0.0001115
wilcox.test(pvals,tooth_sc_pvals)  ##0.3664
wilcox.test(pvals,bone_sc_pvals)   ##0.0001797

wilcox.test(pvals,TF_pvals)       ##0.03356
wilcox.test(pvals,receptors_pvals)##0.144 
wilcox.test(pvals,lig_pvals)      ##0.0002019
wilcox.test(pvals,negReg_pvals)   ##0.06166

pvals <- BC_pvals
diffs <- BC_diffs
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


BC_unique_names <- biteCode_names[!(biteCode_names %in% HFSC_names)]
HFSC_unique_names <- HFSC_names[!(HFSC_names %in% biteCode_names)]
bc_diffs <- sorted_diffs[sorted_names %in% BC_unique_names]
bc_pvals <- sorted_pvals[sorted_names %in% BC_unique_names]
sc_diffs <- sorted_diffs[sorted_names %in% HFSC_unique_names]
sc_pvals <- sorted_pvals[sorted_names %in% HFSC_unique_names]
tetog_diffs <- sorted_diffs[sorted_names %in% lc_tetog_names]
tetog_pvals <- sorted_pvals[sorted_names %in% lc_tetog_names]
lig_diffs <- sorted_diffs[sorted_names %in% lig_names]
lig_pvals <- sorted_pvals[sorted_names %in% lig_names]
receptors_diffs <- sorted_diffs[sorted_names %in% receptors_names]
receptors_pvals <- sorted_pvals[sorted_names %in% receptors_names]
TF_diffs <- sorted_diffs[sorted_names %in% TF_names]
TF_pvals <- sorted_pvals[sorted_names %in% TF_names]
negReg_diffs <- sorted_diffs[sorted_names %in% negReg_names]
negReg_pvals <- sorted_pvals[sorted_names %in% negReg_names]
bone_sc_diffs <- sorted_diffs[sorted_names %in% bone_sc_names]
bone_sc_pvals <- sorted_pvals[sorted_names %in% bone_sc_names]
tooth_sc_diffs <- sorted_diffs[sorted_names %in% tooth_sc_names]
tooth_sc_pvals <- sorted_pvals[sorted_names %in% tooth_sc_names]
sorted_names[sorted_names %in% TF_names]
sorted_names[sorted_names %in% receptors_names]
sorted_names[sorted_names %in% lig_names]
sorted_names[sorted_names %in% negReg_names]
vioplot(pvals,bc_pvals,tetog_pvals,sc_pvals,drawRect = TRUE,
        names = c('all','BiteCode','TeTog','HFSC'),col='red')
dev.print(pdf, paste(Dropbox,'Miller/figures/Multipop_intron4/BC_tetog_HFSC_pvals_vio.pdf',sep=""))
vioplot(pvals,tetog_pvals,drawRect = TRUE,
        names = c('all','TeTog'),col='red')
dev.print(pdf, paste(Dropbox,'Miller/figures/Multipop_intron4/BC_tetog_pvals_vio.pdf',sep=""))

vioplot(pvals,TF_pvals,receptors_pvals,lig_pvals,negReg_pvals,drawRect = TRUE,
        names = c('all','TFs','Receptors','Ligand','NegReg'),col='red')
dev.print(pdf, paste(Dropbox,'Miller/figures/Multipop_intron4/BC_TF_ligands_pvals_vio.pdf',sep=""))


wilcox.test(pvals,bc_pvals)        ##1.194e-15
wilcox.test(pvals,tetog_pvals)     ##7.027e-13
wilcox.test(pvals,sc_pvals)        ##0.018
wilcox.test(pvals,tooth_sc_pvals)  ##0.6202
wilcox.test(pvals,bone_sc_pvals)   ##0.1881

wilcox.test(pvals,TF_pvals)       ##0.06247
wilcox.test(pvals,receptors_pvals)##0.007145
wilcox.test(pvals,lig_pvals)      ##4.668e-07
wilcox.test(pvals,negReg_pvals)   ##0.0004016

wilcox.test(TF_pvals,receptors_pvals) ##0.0846
wilcox.test(TF_pvals,lig_pvals)       ##0.0001402
wilcox.test(TF_pvals,negReg_pvals)    ##0.001727

wilcox.test(receptors_pvals,lig_pvals)    ##0.05291
wilcox.test(receptors_pvals,negReg_pvals) ##0.01607

wilcox.test(lig_pvals,negReg_pvals) ##0.1284

bc_exp <- select_norm_exp(c(biteCode_names),RNA_seq_exp_filtered)
filtered_cor_gset <- filtered_cor[c(colnames(bc_exp),'bmp6','PITX2 (2 of 2)','plod2'),c(colnames(bc_exp),'bmp6','PITX2 (2 of 2)','plod2')]
diag(filtered_cor_gset) <- 0
vioplot(
  filtered_cor[,'bmp6'],filtered_cor_gset[,'bmp6'],
  filtered_cor[,'PITX2 (2 of 2)'],filtered_cor_gset[,'PITX2 (2 of 2)'],
  filtered_cor[,'plod2'],filtered_cor_gset[,'plod2'],
  drawRect = TRUE,names=c('Bmp6 All','Bmp6 BC',
                          'Pitx2 All','Pitx2 BC',
                          'Plod2 All','Plod2 BC'),col=c('green'))


beanplot(filtered_cor[,'bmp6'],filtered_cor_gset[,'bmp6'],
         filtered_cor[,'PITX2 (2 of 2)'],filtered_cor_gset[,'PITX2 (2 of 2)'],
         filtered_cor[,'plod2'],filtered_cor_gset[,'plod2'],
         side='both',names = c('Bmp6','Pitx2','Plod2'),col = list("lightblue", c("red", "black")),ll=NA)
dev.print(pdf, paste(Dropbox,'Miller/figures/Multipop_intron4/Bmp6_Pitx2_Plod2_splitvio.pdf',sep=""))

wilcox.test( filtered_cor[,'bmp6'],filtered_cor_gset[,'bmp6'])$p.value
wilcox.test( filtered_cor[,'PITX2 (2 of 2)'],filtered_cor_gset[,'PITX2 (2 of 2)'])$p.value
wilcox.test( filtered_cor[,'plod2'],filtered_cor_gset[,'plod2'])
