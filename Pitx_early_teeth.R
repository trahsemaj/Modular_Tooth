library(ggplot2)

Dropbox = 'C:/Users/trahs/Dropbox/'
Dropbox = '/home/james/Dropbox/'
paste(Dropbox,'',sep="")


pitx2_data = read.csv(paste(Dropbox,'Miller/data/Mutants/Pitx2/Pitx2_early_teeth.csv',sep=""))
pitx2_data <- pitx2_data[!is.na(pitx2_data$GT),]
pitx2_data <- pitx2_data[!is.na(pitx2_data$TL),]

pitx2_data$GT <- factor(pitx2_data$GT,levels=c('DD','WD','WW'))
summary(lm(TVTP ~ GT + TL + Sex,data=pitx2_data))
summary(lm(TOT ~ GT + TL + Sex,data=pitx2_data))
summary(lm(Dentary ~ GT + TL + Sex,data=pitx2_data))
summary(lm(Premax ~ GT + TL + Sex,data=pitx2_data))

pitx2_data$GT <- factor(pitx2_data$GT,levels=c('WW','WD','DD'))
tooth_lm<-lm(TVTP ~ GT+TL,data=pitx2_data,na.action=na.exclude)
tooth_predict <- predict(tooth_lm,interval="confidence",na.action=na.exclude)
colnames(tooth_predict) <- c('tvtp_fit','tvtp_lwr','tvtp_upr')
pitx2_data<-cbind(pitx2_data,tooth_predict)
ggplot(pitx2_data[!is.na(pitx2_data$TVTP),], aes(x=TL, y=TVTP,color=GT)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(x=TL,y=tvtp_fit))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/Pitx2_early_TVTP_sl.pdf',sep=""))
pitx2_data$TVTP_resid <- resid(lm(TVTP ~ TL,data=pitx2_data,na.action=na.exclude))
ggplot(pitx2_data,aes(x=GT,y=TVTP_resid,color=GT)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Tooth Number') +
  theme(legend.text=element_text(size=8)) + 
  theme(text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Pitx2 \nGenotype",override.aes = list(size=5,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/Pitx2_early_TVTP_box_resid.pdf',sep=""))


AIC(lm(Premax ~ GT+TL,data=pitx2_data))
AIC(lm(Premax ~ GT*TL,data=pitx2_data))


tooth_lm<-lm(Premax ~ GT+TL,data=pitx2_data,na.action=na.exclude)
tooth_predict <- predict(tooth_lm,interval="confidence",na.action=na.exclude)
colnames(tooth_predict) <- c('pm_fit','pm_lwr','pm_upr')
pitx2_data<-cbind(pitx2_data,tooth_predict)
ggplot(pitx2_data, aes(x=TL, y=Premax,color=GT)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))  +
  geom_line(aes(x=TL,y=pm_fit))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/Pitx2_early_pm_sl.pdf',sep=""))


tooth_lm<-lm(Dentary ~ GT+TL,data=pitx2_data,na.action=na.exclude)
tooth_predict <- predict(tooth_lm,interval="confidence",na.action=na.exclude)
colnames(tooth_predict) <- c('den_fit','den_lwr','den_upr')
pitx2_data<-cbind(pitx2_data,tooth_predict)
ggplot(pitx2_data[!is.na(pitx2_data$Dentary),], aes(x=TL, y=Dentary,color=GT)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(x=TL,y=den_fit))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/Pitx2_early_den_sl.pdf',sep=""))
summary(lm(TVTP~TL+GT,data=pitx2_data))
summary(lm(Dentary~TL+GT,data=pitx2_data))
summary(lm(Premax~TL+GT,data=pitx2_data))

summary(lm(Dentary~TL+GT+DPF,data=pitx2_data))
summary(lm(TL~GT,data=pitx2_data))

pitx2_gt_data = read.csv(paste(Dropbox,'Miller/data/Mutants/Pitx2/Pitx2_early_GTs_comb.csv',sep=""))
pitx2_gt_data$total <- pitx2_gt_data$WW + pitx2_gt_data$W1 + pitx2_gt_data$X11 
pitx2_gt_data$obs11 <- pitx2_gt_data$X11/pitx2_gt_data$total

ggplot(pitx2_gt_data, aes(x=Dpf,y=obs11)) + 
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple', "11" = "red")) +
  geom_point(size=3) + 
  geom_line(size=1) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))


