library(ggplot2)

Dropbox = 'C:/Users/trahs/Dropbox/'
Dropbox = '/home/james/Dropbox/'
paste(Dropbox,'',sep="")

plod2_data = read.csv(paste(Dropbox,'Miller/data/Mutants/Plod2/Plod2-comb.csv',sep=""))
##one WT allele doesn't amplify
plod2_data[which(plod2_data$GT == 'HH'),]$GT <- 'HL'
plod2_data[which(plod2_data$GT == 'HL'),]$GT <- 'WD'
plod2_data[which(plod2_data$GT == 'LL'),]$GT <- 'WW'
plod2_data$GT <- factor(plod2_data$GT,levels=c('WW','WD','DD'))
#plod2_data$GT <- factor(plod2_data$GT,levels=c('DD','WD','WW'))

#filter NAs and outliers
plod2_data <- plod2_data[!is.na(plod2_data$GT),]
plod2_data <- plod2_data[!is.na(plod2_data$TVTP),]
plod2_data <- plod2_data[plod2_data$SL>15,]

plod2_data$Clutch <- factor(plod2_data$Clutch)
AIC(lm(TVTP ~ GT+SL,data=plod2_data))
AIC(lm(TVTP ~ GT + Clutch,data=plod2_data))
AIC(lm(TVTP ~ GT+SL+Clutch,data=plod2_data))
##This still seems to be the best model
AIC(lm(TVTP ~ GT+SL:Clutch,data=plod2_data))
AIC(lm(TVTP ~ GT+SL*Clutch,data=plod2_data))


summary(  lm(TVTP ~ GT+SL:Clutch,data=plod2_data))

summary(  lm(SL ~ GT+Clutch,data=plod2_data))





tooth_lm<-lm(TVTP ~ GT+SL:Clutch,data=plod2_data)
tooth_predict <- predict(tooth_lm,interval="confidence")
colnames(tooth_predict) <- c('TVTP_fit','TVTP_low','TVTP_up')
plod2_data<-cbind(plod2_data,tooth_predict)
ggplot(plod2_data, aes(x=SL, y=TVTP,color=GT,linetype=Clutch,shape=Clutch)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=TVTP_fit),size=2)

ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_allGT_TVTP_SL.pdf',sep=""))

plod2_data$tvtp_res <- resid(lm(TVTP ~ SL:Clutch,data=plod2_data))
qplot(GT,tvtp_res, geom=c("boxplot","jitter"),color=GT,data = plod2_data)  +
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  labs(x='GT',y='Corrected Tooth Number',color='Plod2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 


ggplot(plod2_data,aes(x=GT,y=tvtp_res,color=GT)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Tooth Number') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Plod2 \nGenotype",override.aes = list(size=3,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_allGT_tvtp_res_box.pdf',sep=""))

summary(  lm(SL ~ GT+Clutch,data=plod2_data))

sl_lm <- lm(SL ~ Clutch,data=plod2_data)
plod2_data$sl_res <- resid(sl_lm)

qplot(GT,sl_res, geom=c("boxplot","jitter"),color=GT,data = plod2_data)  +
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  labs(x='GT',y='Corrected_Standard Length',color='Plod2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 
summary(lm(SL~GT+Clutch,data=plod2_data))


ggplot(plod2_data[plod2_data$sl_res>-10,],aes(x=GT,y=sl_res,color=GT)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Length') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Plod2 \nGenotype",override.aes = list(size=3,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_allGT_sl_res_box.pdf',sep=""))

################
##Break it down
##################
with(plod2_3gt, table(GT, Clutch))
plod2_3gt <- plod2_data[!(plod2_data$Clutch %in% c('TG1346A','TG1346B')),]
#plod2_3gt <- plod2_3gt[!is.na(plod2_3gt$Sex),]

AIC(lm(TVTP ~ GT+SL,data=plod2_3gt))
AIC(lm(TVTP ~ GT+Clutch,data=plod2_3gt))
AIC(lm(TVTP ~ GT+SL+Clutch,data=plod2_3gt))
AIC(lm(TVTP ~ GT+SL:Clutch,data=plod2_3gt))
AIC(lm(TVTP ~ GT+SL*Clutch,data=plod2_3gt))
AIC(lm(TVTP ~ GT:SL+SL:Clutch,data=plod2_3gt))
AIC(lm(TVTP ~ GT*SL*Clutch,data=plod2_3gt))

plod2_3gt$GT <- factor(plod2_3gt$GT,levels=c('DD','WD','WW'))
summary(lm(TVTP ~ GT+SL:Clutch,data=plod2_3gt))
plod2_3gt$GT <- factor(plod2_3gt$GT,levels=c('WW','WD','DD'))
summary(lm(TVTP ~ GT+Clutch,data=plod2_3gt))



tooth_lm<-lm(TVTP ~ GT+SL*Clutch,data=plod2_3gt)
tooth_predict <- predict(tooth_lm,interval="confidence")
colnames(tooth_predict) <- c('tvtp_fit','tvtp_lwr','tvtp_upr')
plod2_3gt<-cbind(plod2_3gt,tooth_predict)
ggplot(plod2_3gt, aes(x=SL, y=TVTP,color=GT,linetype=Clutch,shape=Clutch)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=tvtp_fit),size=2)

ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_TVTP_SL.pdf',sep=""))

plod2_3gt$tvtp_res <- resid(lm(TVTP ~ SL:Clutch,data=plod2_3gt,na.action=na.exclude))
ggplot(plod2_3gt,aes(x=GT,y=tvtp_res,color=GT)) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Length',x='Plod2 Genotype') +
  theme(legend.text=element_text(size=8)) + 
  theme(text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Plod2 \nGenotype",override.aes = list(size=5,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_tvtp_res_box.pdf',sep=""))


AIC(lm(TDTP1 ~ GT+SL,data=plod2_3gt))
AIC(lm(TDTP1 ~ GT+SL+Sex,data=plod2_3gt))

summary(  lm(TDTP1 ~ GT+SL,data=plod2_3gt))

tooth_lm<-lm(TDTP1 ~ GT+SL,data=plod2_3gt,na.action=na.exclude)
tooth_predict <- predict(tooth_lm,interval="confidence",na.action=na.exclude)
colnames(tooth_predict) <- c('tdtp1_fit','tdtp1_lwr','tdtp1_upr')
plod2_3gt<-cbind(plod2_3gt,tooth_predict)
ggplot(plod2_3gt[!is.na(plod2_3gt$TDTP1),], aes(x=SL, y=TDTP1,color=GT)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  xlim(c(15,21)) +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=tdtp1_fit),size=2)

ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_TDTP1_SL.pdf',sep=""))

plod2_3gt$tdtp1_res <- resid(lm(TDTP1 ~ SL,data=plod2_3gt,na.action=na.exclude))
ggplot(plod2_3gt,aes(x=GT,y=tdtp1_res,color=GT)) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Tooth Number',x='Plod2 Genotype') +
  theme(legend.text=element_text(size=8)) + 
  theme(text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Plod2 \nGenotype",override.aes = list(size=5,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_tdtp1_res_box.pdf',sep=""))

AIC(lm(TDTP2 ~ GT+SL,data=plod2_3gt))
AIC(lm(TDTP2 ~ GT+SL+Sex,data=plod2_3gt))
AIC(lm(TDTP2 ~ GT+SL:Sex,data=plod2_3gt))
AIC(lm(TDTP2 ~ GT+SL*Sex,data=plod2_3gt))
summary(  lm(TDTP2 ~ GT+SL+Sex,data=plod2_3gt))

tooth_lm<-lm(TDTP2 ~ GT+SL+Sex,data=plod2_3gt,na.action=na.exclude)
tooth_predict <- predict(tooth_lm,interval="confidence",na.action=na.exclude)
colnames(tooth_predict) <- c('tdtp2_fit','tdtp2_lwr','tdtp2_upr')
plod2_3gt<-cbind(plod2_3gt,tooth_predict)
ggplot(plod2_3gt[!is.na(plod2_3gt$TDTP2),], aes(x=SL, y=TDTP2,color=GT,shape=Sex)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  xlim(c(15,21)) +
  geom_line(aes(x=SL,y=tdtp2_fit,color=GT,linetype=Sex),size=1,inherit.aes=FALSE)

ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_TDTP2_SL.pdf',sep=""))

plod2_3gt$tdtp2_res <- resid(lm(TDTP2 ~ SL+Sex,data=plod2_3gt,na.action=na.exclude))
ggplot(plod2_3gt,aes(x=GT,y=tdtp2_res,color=GT)) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Tooth Number',x='Plod2 Genotype') +
  theme(legend.text=element_text(size=8)) + 
  theme(text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Plod2 \nGenotype",override.aes = list(size=5,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_tdtp2_res_box.pdf',sep=""))

AIC(lm(TotalDentary ~ GT+SL,data=plod2_3gt))
AIC(lm(TotalDentary ~ GT+SL+Sex,data=plod2_3gt))
AIC(lm(TotalDentary ~ GT+SL:Sex,data=plod2_3gt))
AIC(lm(TotalDentary ~ GT+SL*Sex,data=plod2_3gt))
summary(  lm(TotalDentary ~ GT+SL+Sex,data=plod2_3gt))

tooth_lm<-lm(TotalDentary ~ GT+SL+Sex,data=plod2_3gt,na.action=na.exclude)
tooth_predict <- predict(tooth_lm,interval="confidence",na.action=na.exclude)
colnames(tooth_predict) <- c('den_fit','den_lwr','den_upr')
plod2_3gt<-cbind(plod2_3gt,tooth_predict)
ggplot(plod2_3gt, aes(x=SL, y=TotalDentary,color=GT,shape=Sex)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  xlim(c(15,21)) +
  geom_line(aes(x=SL,y=den_fit,color=GT,linetype=Sex),size=1,inherit.aes=FALSE)

ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_Dentary_SL.pdf',sep=""))

plod2_3gt$den_res <- resid(lm(TotalDentary ~ SL+Sex,data=plod2_3gt,na.action=na.exclude))
ggplot(plod2_3gt,aes(x=GT,y=den_res,color=GT)) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Tooth Number',x='Plod2 Genotype') +
  theme(legend.text=element_text(size=8)) + 
  theme(text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Plod2 \nGenotype",override.aes = list(size=5,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_dentary_res_box.pdf',sep=""))

plod2_3gt$GT <- factor(plod2_3gt$GT,levels=c('DD','WD','WW'))
summary(  lm(TotalPreMax ~ GT+SL+Sex,data=plod2_3gt,na.action=na.exclude))
plod2_3gt$GT <- factor(plod2_3gt$GT,levels=c('WW','WD','DD'))


AIC(lm(TotalPreMax ~ GT+SL,data=plod2_3gt))
AIC(lm(TotalPreMax ~ GT+SL+Sex,data=plod2_3gt))
AIC(lm(TotalPreMax ~ GT+SL:Sex,data=plod2_3gt))


tooth_lm<-lm(TotalPreMax ~ GT+SL+Sex,data=plod2_3gt,na.action=na.exclude)
tooth_predict <- predict(tooth_lm,interval="confidence",na.action=na.exclude)
colnames(tooth_predict) <- c('pm_fit','pm_lwr','pm_upr')
plod2_3gt<-cbind(plod2_3gt,tooth_predict)
ggplot(plod2_3gt, aes(x=SL, y=TotalPreMax,color=GT,shape=Sex)) + 
  scale_colour_manual(values =c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  xlim(c(15,21)) +
  geom_line(aes(x=SL,y=pm_fit,color=GT,linetype=Sex),size=1,inherit.aes=FALSE)

ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_Premax_SL.pdf',sep=""))

plod2_3gt$pm_res <- resid(lm(TotalPreMax ~ SL+Sex,data=plod2_3gt,na.action=na.exclude))
ggplot(plod2_3gt,aes(x=GT,y=pm_res,color=GT)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Tooth Number',x='Plod2 Genotype') +
  theme(legend.text=element_text(size=8)) + 
  theme(text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Plod2 \nGenotype",override.aes = list(size=5,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_Premax_res_box.pdf',sep=""))

plod2_3gt$GT <- factor(plod2_3gt$GT,levels=c('DD','WD','WW'))
summary(lm(SL ~ GT+Clutch,data=plod2_3gt,na.action=na.exclude))
plod2_3gt$GT <- factor(plod2_3gt$GT,levels=c('WW','WD','DD'))
plod2_3gt$sl_res <- resid(lm(SL ~ Clutch,data=plod2_3gt,na.action=na.exclude))
ggplot(plod2_3gt,aes(x=GT,y=sl_res,color=GT)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Length',x='Plod2 Genotype') +
  theme(legend.text=element_text(size=8)) + 
  theme(text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Plod2 \nGenotype",override.aes = list(size=5,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_SL_res_box.pdf',sep=""))

qplot(GT,SL, geom=c("boxplot","jitter"),color=GT,data = plod2_3gt)  +
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  labs(x='GT',y='Standard Length',color='Plod2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 
summary(lm(SL~GT+Clutch,data=plod2_data))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Plod2/Plod2_3GT_sl_box.pdf',sep=""))
