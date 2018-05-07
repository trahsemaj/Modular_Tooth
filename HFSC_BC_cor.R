library(ggplot2)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(FactoMineR)

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

  #select_exp <- exp[toupper(rownames(exp)) %in% nameList,]
  select_exp <-exp[grep(paste(nameList,collapse='|'),rownames(exp),ignore.case=TRUE),]
  #select_exp <- exp[nameList,]
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


##Load HFSC exp
HFSC <- read.table('/home/james/Dropbox/Miller/data/BED/HFSCGS_snames.l',sep='\t')
HFSC_names = as.vector(HFSC[,c(1)])
HFSC_exp <- select_norm_exp(c(HFSC_names),RNA_seq_exp_filtered)
HFSC_cor <- cor(HFSC_exp,method='pearson')
diag(HFSC_cor) <- 0

##Load BC exp
biteCode <- read.table('/home/james/Dropbox/Miller/data/BED/BITECODE.l',sep='\t')
biteCode_names = as.vector(biteCode[,c(1)])
BC_exp <- select_norm_exp(c(biteCode_names),RNA_seq_exp_filtered)
BC_cor <- cor(BC_exp,method='pearson')
diag(BC_cor) <- 0

tetog <- read.table('/home/james/Dropbox/Miller/data/BED/TeToG_snames.l',sep='\t')
tetog_names = as.vector(tetog[,c(1)])
tetog_exp <- select_norm_exp_upper(c(tetog_names),RNA_seq_exp_filtered)
tetog_cor <- cor(tetog_exp,method='pearson')
lc_tetog_names <- colnames(tetog_exp)
lc_tetog_names

tetog_epi <- read.table('/home/james/Dropbox/Miller/data/BED/TeToG_epi.l',sep='\t')
tetog_epi_names = as.vector(tetog_epi[,c(1)])
tetog_epi_exp <- select_norm_exp_upper(c(tetog_epi_names),RNA_seq_exp_filtered)
tetog_epi_cor <- cor(tetog_epi_exp,method='pearson')
lc_tetog_epi_names <- colnames(tetog_epi_exp)
lc_tetog_epi_names

tetog_mes <- read.table('/home/james/Dropbox/Miller/data/BED/TeToG_mes.l',sep='\t')
tetog_mes_names = as.vector(tetog_mes[,c(1)])
tetog_mes_exp <- select_norm_exp_upper(c(tetog_mes_names),RNA_seq_exp_filtered)
tetog_mes_cor <- cor(tetog_mes_exp,method='pearson')
lc_tetog_mes_names <- colnames(tetog_mes_exp)
lc_tetog_mes_names


plot_bc_hfsc <- function(name){
  
  
  
  bc_pval <- signif(wilcox.test( filtered_cor[,name],BC_cor[,name])$p.value,digits=4)
  hfsc_pval <- signif(wilcox.test( filtered_cor[,name],HFSC_cor[,name])$p.value,digits=4)
  diff_pval <- signif(wilcox.test( BC_cor[,name],HFSC_cor[,name])$p.value,digits=4)
  vioplot(
    filtered_cor[,name],BC_cor[,name],HFSC_cor[,name],
    drawRect = TRUE,names=c(paste('All',name,diff_pval,sep = '_'),
                            paste('BC',bc_pval,sep = '_'),
                            paste('HFSC',hfsc_pval,sep = '_'))
  )
  return(c(wilcox.test( filtered_cor[,name],BC_cor[,name])$p.value,wilcox.test( filtered_cor[,name],HFSC_cor[,name])$p.value))
}


plot_tetog_epi_mes <- function(name){
  
  tetog_exp <- select_norm_exp(c(lc_tetog_names,name),RNA_seq_exp_filtered)
  tetog_cor <- cor(tetog_exp,method='pearson')
  tetog_epi_exp <- select_norm_exp(c(lc_tetog_epi_names,name),RNA_seq_exp_filtered)
  tetog_epi_cor <- cor(tetog_epi_exp,method='pearson')
  tetog_mes_exp <- select_norm_exp(c(lc_tetog_mes_names,name),RNA_seq_exp_filtered)
  tetog_mes_cor <- cor(tetog_mes_exp,method='pearson')
  
  tetog_pval <- signif(wilcox.test( filtered_cor[,name],tetog_cor[,name])$p.value,digits=4)
  epi_pval <- signif(wilcox.test( filtered_cor[,name],tetog_epi_cor[,name])$p.value,digits=4)
  mes_pval <- signif(wilcox.test( filtered_cor[,name],tetog_mes_cor[,name])$p.value,digits=4)
  diff_pval <- signif(wilcox.test( tetog_epi_cor[,name],tetog_mes_cor[,name])$p.value,digits=4)
  vioplot(
    filtered_cor[,name],tetog_cor[,name],tetog_epi_cor[,name],tetog_mes_cor[,name],
    drawRect = TRUE,names=c(paste('All',name,diff_pval,sep = '_'),
                            paste('Tetog',tetog_pval,sep = '_'),
                            paste('epi',epi_pval,sep = '_'),
                            paste('mes',mes_pval,sep = '_'))
  )
  return(c(wilcox.test( filtered_cor[,name],tetog_epi_cor[,name])$p.value,wilcox.test( filtered_cor[,name],tetog_mes_cor[,name])$p.value))
}



RNA_seq_exp_filtered[grep('PITX2 (2 of 2)',rownames(RNA_seq_exp_filtered)),]
head(sort(filtered_cor['bmp4',],decreasing=TRUE),n=200)
head(sort(filtered_cor[,'msx2b (2 of 2)'],decreasing=TRUE),n=50)
#lc_tetog_names <- c(lc_tetog_names,'LGR6')
#lc_tetog_names <- c('plod2','PITX2 (2 of 2)','bmp6','ZNF106 (2 of 2)','grem2a','LGR6','LGR4','MSX2','msx2b (1 of 2)','msx2b (2 of 2)')
tested_pvals <- c()
for (n in lc_tetog_names){
  BC_exp <- select_norm_exp(c(biteCode_names,n),RNA_seq_exp_filtered)
  BC_cor <- cor(BC_exp,method='pearson')
  diag(BC_cor) <- 0
  HFSC_exp <- select_norm_exp(c(HFSC_names,n),RNA_seq_exp_filtered)
  HFSC_cor <- cor(HFSC_exp,method='pearson')
  diag(HFSC_cor) <- 0
  pvals <- plot_bc_hfsc(n)
  tested_pvals <- c(tested_pvals,c(n,pvals))
}
tested_pvals




head(sort(filtered_cor['bmp4',],decreasing=TRUE),n=200)
head(sort(filtered_cor[,'msx2b (2 of 2)'],decreasing=TRUE),n=50)
#lc_tetog_names <- c(lc_tetog_names,'LGR6')
names_to_test <- c('plod2','PITX2 (2 of 2)','bmp6','ZNF106 (2 of 2)','grem2a','LGR6','LGR4','MSX2','msx2b (1 of 2)','msx2b (2 of 2)')
tested_pvals <- c()
names_to_test <- lc_tetog_epi_names
for (n in names_to_test){
  pvals <- plot_tetog_epi_mes(n)
  tested_pvals <- c(tested_pvals,c(n,pvals))
}


vioplot(
  filtered_cor[,'bmp6'],BC_cor[,'bmp6'],HFSC_cor[,'bmp6']
)
wilcox.test( filtered_cor[,'bmp6'],BC_cor[,'bmp6'])$p.value
wilcox.test( filtered_cor[,'bmp6'],HFSC_cor[,'bmp6'])$p.value

vioplot(
  filtered_cor[,'LGR6'],BC_cor[,'LGR6'],HFSC_cor[,'LGR6']
)
wilcox.test( filtered_cor[,'LGR6'],BC_cor[,'LGR6'])$p.value
wilcox.test( filtered_cor[,'LGR6'],HFSC_cor[,'LGR6'])$p.value

filtered_cor[,'PITX2 (2 of 2)']
vioplot(
  filtered_cor[,'PITX2 (2 of 2)'],BC_cor[,'PITX2 (2 of 2)'],HFSC_cor[,'PITX2 (2 of 2)']
)
wilcox.test( filtered_cor[,'PITX2 (2 of 2)'],BC_cor[,'PITX2 (2 of 2)'])$p.value
wilcox.test( filtered_cor[,'PITX2 (2 of 2)'],HFSC_cor[,'PITX2 (2 of 2)'])$p.value

vioplot(
  filtered_cor[,'plod2'],BC_cor[,'plod2'],HFSC_cor[,'plod2']
)
wilcox.test( filtered_cor[,'plod2'],BC_cor[,'plod2'])$p.value
wilcox.test( filtered_cor[,'plod2'],HFSC_cor[,'plod2'])$p.value


vioplot(
  filtered_cor[,'shha'],BC_cor[,'shha'],HFSC_cor[,'shha']
)
wilcox.test( filtered_cor[,'shha'],BC_cor[,'shha'])$p.value
wilcox.test( filtered_cor[,'shha'],HFSC_cor[,'shha'])$p.value
