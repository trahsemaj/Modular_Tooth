library(ggplot2)
#install.packages('beeswarm')
library(beeswarm)
all_qtl_data <- read.csv('/home/james/Dropbox/Miller/data/QTL_Data/combined_raw_cross_sheet_080217.csv')
all_qtl_data <- read.csv('C:/Users/trahs/Dropbox/Miller/data/QTL_Data/combined_raw_cross_sheet_080217.csv')
all_qtl_data[which(all_qtl_data$TVTP == 0),'TVTP'] <- NA



#tvtp = all_qtl_data$LVTP + all_qtl_data$RVTP
#all_qtl_data['TVTP'] <- tvtp
#all_qtl_data[is.na(all_qtl_data$LVTP),'TVTP'] <- all_qtl_data[is.na(all_qtl_data$LVTP),]$RVTP * 2
#all_qtl_data[is.na(all_qtl_data$RVTP),'TVTP'] <- all_qtl_data[is.na(all_qtl_data$RVTP),]$LVTP * 2

#all_qtl_data <- all_qtl_data[!is.na(all_qtl_data$BMP6_raw),]
plot(TVTP ~ BMP6_raw, data=all_qtl_data)
all_qtl_data['Bmp6_MF'] <- all_qtl_data$BMP6_raw
all_qtl_data[which(all_qtl_data$Bmp6_MF == 'notMM'),'Bmp6_MF'] <- NA
all_qtl_data[which(all_qtl_data$Bmp6_MF == 'notFF'),'Bmp6_MF'] <- NA

CCD_data <- all_qtl_data[which(all_qtl_data$Cross == 'CCDxLITC'),]
table(CCD_data$BMP6_raw)
TVTP_SL_res <- resid(lm(TVTP ~ SL,data=CCD_data,na.action=na.exclude))
CCD_data['TVTP_SL_res'] <- TVTP_SL_res
plot(TVTP_SL_res ~BMP6_raw,data=CCD_data)
summary(aov(TVTP_SL_res ~BMP6_raw,data=CCD_data))


plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'BEPAxLITC'),])
all_qtl_data[which(all_qtl_data$Cross == 'BEPAxLITC' &  all_qtl_data$BMP6_raw == 'A'),]$Bmp6_MF <- 'MM'
all_qtl_data[which(all_qtl_data$Cross == 'BEPAxLITC' &  all_qtl_data$BMP6_raw == 'B'),]$Bmp6_MF <- 'FF'
all_qtl_data[which(all_qtl_data$Cross == 'BEPAxLITC' &  all_qtl_data$BMP6_raw == 'H'),]$Bmp6_MF <- 'MF'

##B appears to be F, should be additive?
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'FTCxLITC'),])
all_qtl_data[which(all_qtl_data$Cross == 'FTCxLITC' &  all_qtl_data$BMP6_raw == 'A'),]$Bmp6_MF <- 'MM'
all_qtl_data[which(all_qtl_data$Cross == 'FTCxLITC' &  all_qtl_data$BMP6_raw == 'B'),]$Bmp6_MF <- 'FF'
all_qtl_data[which(all_qtl_data$Cross == 'FTCxLITC' &  all_qtl_data$BMP6_raw == 'H'),]$Bmp6_MF <- 'MF'

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'CERCxLITC'),])

##D RFLP
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'CCDxLITC'),])
all_qtl_data[which(all_qtl_data$Cross == 'CCDxLITC' &  all_qtl_data$BMP6_raw == 'DD'),]$Bmp6_MF <- 'FF'
all_qtl_data[which(all_qtl_data$Cross == 'CCDxLITC' &  all_qtl_data$BMP6_raw == 'MD'),]$Bmp6_MF <- 'MF'

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'LITCxPAXB'),])


##chrXXI H = CERC, L = Marine
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'CERCxRABS'),])
all_qtl_data[which(all_qtl_data$Cross == 'CERCxRABS' &  all_qtl_data$BMP6_raw == 'HH'),]$Bmp6_MF <- 'MM'
all_qtl_data[which(all_qtl_data$Cross == 'CERCxRABS' &  all_qtl_data$BMP6_raw == 'HL'),]$Bmp6_MF <- 'FF'
all_qtl_data[which(all_qtl_data$Cross == 'CERCxRABS' &  all_qtl_data$BMP6_raw == 'LL'),]$Bmp6_MF <- 'MF'

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_L'),])

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_A'),])

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxRABS'),])


##Set b and d to F, as BD = DD, even though MB looks different than MD
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_WV'),])
all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_WV' &  all_qtl_data$BMP6_raw == 'bb'),]$Bmp6_MF <- 'FF'
all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_WV' &  all_qtl_data$BMP6_raw == 'bd'),]$Bmp6_MF <- 'FF'
all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_WV' &  all_qtl_data$BMP6_raw == 'mb'),]$Bmp6_MF <- 'MF'
all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_WV' &  all_qtl_data$BMP6_raw == 'md'),]$Bmp6_MF <- 'MF'
all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_WV' &  all_qtl_data$BMP6_raw == 'mm'),]$Bmp6_MF <- 'MM'


##Set d to F and b to M
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxLITC_28'),])
points(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxLITC_28'),])
all_qtl_data[which(all_qtl_data$Cross == 'PAXBxLITC_28' &  all_qtl_data$BMP6_raw == 'dd'),]$Bmp6_MF <- 'FF'
all_qtl_data[which(all_qtl_data$Cross == 'PAXBxLITC_28' &  all_qtl_data$BMP6_raw == 'Md'),]$Bmp6_MF <- 'MF'
all_qtl_data[which(all_qtl_data$Cross == 'PAXBxLITC_28' &  all_qtl_data$BMP6_raw == 'Mb'),]$Bmp6_MF <- 'MM'




##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxLITC_29'),])

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'LITCxENOS'),])

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'LITCxPRIB'),])

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'HUTUxLITC'),])

##all set
plot(TVTP ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'LITCxFTC'),])


all_qtl_data$Bmp6_MF <- factor(all_qtl_data$Bmp6_MF)
#all_qtl_data <- all_qtl_data[!is.na(all_qtl_data$Bmp6_MF),]

plot(TVTP ~ Bmp6_MF, data=all_qtl_data)
qplot(Bmp6_MF,TVTP, geom=c("boxplot"),color=Cross,data = all_qtl_data)  + theme(text = element_text(size=20))


summary(lm(TVTP ~ SL,data=all_qtl_data))
summary(lm(TVTP ~ SL + Tank,data=all_qtl_data))
summary(lm(TVTP ~ SL + Tank + Cross,data=all_qtl_data))
length(all_qtl_data$SL)

AIC(lm(TVTP ~ SL,data=all_qtl_data))
AIC(lm(TVTP ~ SL + Tank,data=all_qtl_data))
AIC(lm(TVTP ~ SL + Tank + Cross,data=all_qtl_data))
AIC(lm(TVTP ~ SL*Cross + Tank,data=all_qtl_data))
AIC(lm(TVTP ~ SL:Cross + Tank,data=all_qtl_data))
AIC(lm(TVTP ~ SL:Cross + Tank+SL,data=all_qtl_data))
AIC(lm(TVTP ~ SL:Tank + Tank+SL,data=all_qtl_data))
AIC(lm(TVTP ~ SL*Cross,data=all_qtl_data))
#AIC(lm(TVTP ~ SL*Tank*Cross,data=all_qtl_data))
AIC(lm(TVTP ~ SL+SL:Tank,data=all_qtl_data))


summary(lm(TVTP ~ SL:Cross + Tank + Cross:Bmp6_MF,data=all_qtl_data))

table(CCD_data$Bmp6_MF)
tvtp_res <- resid(lm(TVTP ~ SL:Cross + Tank,data=all_qtl_data,na.action=na.exclude))
all_qtl_data$tvtp_res <- tvtp_res

summary(aov(lm(TVTP ~ SL:Cross + Tank + Bmp6_MF,data=all_qtl_data,na.action=na.exclude)))
qplot(Bmp6_MF,TVTP, geom=c("boxplot"),color=Cross,data = all_qtl_data)  + theme(text = element_text(size=20))

qplot(Bmp6_MF,tvtp_res, geom=c("boxplot"),color=Cross,data = all_qtl_data)  + theme(text = element_text(size=20))



all_qtl_data <- all_qtl_data[!is.na(all_qtl_data$Bmp6_MF),]

mean_tvtp <- aggregate(tvtp_res~Cross, FUN=mean,data=all_qtl_data[which(all_qtl_data$Bmp6_MF == 'MM'),])
mean_tvtp['tvtp_res']

summary(lm(tvtp_res ~ Cross:Bmp6_MF,data=all_qtl_data[which(all_qtl_data$Bmp6_MF != 'MF'),]))

head(sapply(all_qtl_data$Cross,function(x) mean(all_qtl_data[which(all_qtl_data$Cross == x & all_qtl_data$Bmp6_MF == 'MM'),'tvtp_res'])))
##get mean diff MM vs FF
all_qtl_data$Mean_MM <- sapply(all_qtl_data$Cross,function(x) mean(all_qtl_data[which(all_qtl_data$Cross == x & all_qtl_data$Bmp6_MF == 'MM'),'tvtp_res'],na.rm=TRUE))
all_qtl_data$tvtp_res_norm <- all_qtl_data$tvtp_res - all_qtl_data$Mean_MM

##get the pval for the aov for each cross
all_qtl_data$GT_pval <- sapply(all_qtl_data$Cross,function(x) summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == x),]))[[1]]$'Pr(>F)'[[1]])
#summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == x),]))[[1]]$'Pr(>F)'[[1]]




FF_norm <- aggregate(tvtp_res_norm~Cross, FUN=mean,data=all_qtl_data[which(all_qtl_data$Bmp6_MF == 'FF'),])
#levels(all_qtl_data$Cross) <- as.vector(FF_norm[with(FF_norm, order(tvtp_res_norm)),]$Cross)
as.vector(FF_norm[with(FF_norm, order(tvtp_res_norm)),]$Cross)
##Plot order based off norm FF, MM difference
all_qtl_data$plot_order <- factor(all_qtl_data$Cross, as.vector(FF_norm[with(FF_norm, order(tvtp_res_norm)),]$Cross))
plot_order <- as.vector(FF_norm[with(FF_norm, order(tvtp_res_norm)),]$Cross)

pval_arg <- aggregate(GT_pval~Cross, FUN=mean,data=all_qtl_data)
pval_arg
#levels(all_qtl_data$Cross) <- as.vector(FF_norm[with(FF_norm, order(tvtp_res_norm)),]$Cross)
as.vector(pval_arg[with(pval_arg, order(GT_pval)),]$Cross)
##plot order based on pvals
plot_order <- as.vector(pval_arg[with(pval_arg, order(GT_pval)),]$Cross)
all_qtl_data$plot_order <- factor(all_qtl_data$Cross, as.vector(pval_arg[with(pval_arg, order(GT_pval,decreasing=TRUE)),]$Cross))

plot_order
##in LITC29, 29=BB????
alpha=.1
tested_hyp = length(plot_order)
HB_cutoff <- c(1:tested_hyp)
HB_cutoff <- alpha/(tested_hyp + 1 - HB_cutoff)
plot_order[pval_arg[with(pval_arg, order(GT_pval)),]$GT_pval < HB_cutoff]

alpha=.05
tested_hyp = length(plot_order)
HB_cutoff <- c(1:tested_hyp)
HB_cutoff <- alpha/(tested_hyp + 1 - HB_cutoff)
plot_order[pval_arg[with(pval_arg, order(GT_pval)),]$GT_pval < HB_cutoff]


alpha=.01
tested_hyp = length(plot_order)
HB_cutoff <- c(1:tested_hyp)
HB_cutoff <- alpha/(tested_hyp + 1 - HB_cutoff)
plot_order[pval_arg[with(pval_arg, order(GT_pval)),]$GT_pval < HB_cutoff]


alpha=.001
tested_hyp = length(plot_order)
HB_cutoff <- c(1:tested_hyp)
HB_cutoff <- alpha/(tested_hyp + 1 - HB_cutoff)
plot_order[pval_arg[with(pval_arg, order(GT_pval)),]$GT_pval < HB_cutoff]

all_qtl_data$Bmp6_MF <- factor(all_qtl_data$Bmp6_MF,levels=c('MM','MF','FF'))
#write.table(all_qtl_data,file='/home/james/Dropbox/Miller/data/QTL_Data/combined_cross_sheet_forR_032218.csv',sep='\t',quote=FALSE)
##FTCxLITC markers are a bit off??
##HB correction leaves BEPA as barely significant (.0123 < .0125)
qplot(plot_order,tvtp_res_norm, geom=c("boxplot"),color=Bmp6_MF,data = all_qtl_data,outlier.shape = 20)  +
  scale_colour_manual(values = c("FF" = "blue", 'MF' = 'darkgreen', "MM" = "red")) +
  labs(x='Cross',y='Corrected Tooth Number',color='Bmp6 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 

ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/chr21_Tooth_num_boxplot.pdf',width=6,height=4.5,units='in')


qplot(plot_order,tvtp_res_norm, geom=c("boxplot"),color=Bmp6_MF,data = all_qtl_data[which(all_qtl_data$Cross == 'PAXBxRABS'),],outlier.shape = 20)  +
  scale_colour_manual(values = c("FF" = "blue", 'MF' = 'darkgreen', "MM" = "red")) +
  labs(x='Cross',y='Corrected Tooth Number',color='Bmp6 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/chr21_Tooth_num_boxplot_PxR.pdf',width=6,height=4.5,units='in')

qplot(Bmp6_MF,tvtp_res_norm, geom=c("boxplot","jitter"),color=Bmp6_MF,data = all_qtl_data[which(all_qtl_data$Cross == 'PAXBxRABS'),],outlier.shape = 20)  +
  scale_colour_manual(values = c("FF" = "blue", 'MF' = 'darkgreen', "MM" = "red")) +
  labs(x='Cross',y='Corrected Tooth Number',color='Bmp6 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/chr21_Tooth_num_boxplot_jitter_PxR.pdf',width=6,height=4.5,units='in')


summary(aov(tvtp_res_norm ~ Bmp6_MF,data = all_qtl_data[which(all_qtl_data$Cross == 'PAXBxRABS'),]))
TukeyHSD(aov(tvtp_res_norm ~ Bmp6_MF,data = all_qtl_data[which(all_qtl_data$Cross == 'PAXBxRABS'),]))

ggplot(all_qtl_data, aes(x=SL, y=TVTP,color=Cross)) + 
  geom_point(size=.5) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(all_qtl_data, aes(x=SL, y=TVTP,color=Tank)) + 
  geom_point(size=.5) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)
table(all_qtl_data$Cross)
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'BEPAxLITC'),]))[[1]]$'Pr(>F)'[[1]]
beeswarm(tvtp_res_norm~plot_order,data = all_qtl_data,pch=16,pwcol = Bmp6_MF)
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'BEPAxLITC'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'FTCxLITC'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'CCDxLITC'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'LITCxPAXB'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'CERCxRABS'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'CERCxLITC'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_L'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_A'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxJAMA_WV'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxLITC_28'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'PAXBxLITC_29'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'LITCxENOS'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'LITCxPRIB'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'HUTUxLITC'),]))
summary(aov(tvtp_res ~ Bmp6_MF, data=all_qtl_data[which(all_qtl_data$Cross == 'LITCxFTC'),]))
