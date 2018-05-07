library(ggplot2)

Dropbox = 'C:/Users/trahs/Dropbox/'
Dropbox = '/home/james/Dropbox/'
paste(Dropbox,'',sep="")
##PITX2 comb
##NTC Pitx2
#install.packages('ggplot2')
#pitx2_data = read.csv('C:/Users/trahs/Dropbox/Miller/data/Mutants/Pitx2/NTC_PITX2_081817.csv')
pitx2_data = read.csv(paste(Dropbox,'Miller/data/Mutants/Pitx2/NTC_PITX2_081817.csv',sep=""))
factor(pitx2_data$Genotype)
pitx2_data$GT <- factor(pitx2_data$Genotype,levels=c('WW','W37','W2','W1','3737'))

clean_pitx2_data <- pitx2_data[!is.na(pitx2_data$Genotype),]
clean_pitx2_data <- clean_pitx2_data[!is.na(clean_pitx2_data$TVTP),]
clean_pitx2_data <- clean_pitx2_data[!is.na(clean_pitx2_data$SL),]
##no effect on Spine Length
ggplot(pitx2_data, aes(x=SL, y=Pelv_Left,color=GT)) + 
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'green','W2'='green','W1'='purple','3737'='red')) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_allGT_Pelv_Left.pdf',sep=""))
ggplot(pitx2_data, aes(x=SL, y=Pelv_Right,color=GT)) + 
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'green','W2'='green','W1'='purple','3737'='red')) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_allGT_Pelv_Right.pdf',sep=""))



AIC(lm(TVTP ~ GT+SL,data=clean_pitx2_data))
AIC(lm(TVTP ~ GT+Clutch,data=clean_pitx2_data))
AIC(lm(TVTP ~ GT+SL+Clutch,data=clean_pitx2_data))
AIC(lm(TVTP ~ GT+SL:Clutch,data=clean_pitx2_data))
AIC(lm(TVTP ~ GT+SL*Clutch,data=clean_pitx2_data))

summary(lm(TVTP ~ GT+SL,data=clean_pitx2_data))
summary(lm(TVTP ~ GT+Clutch,data=clean_pitx2_data))
summary(lm(TVTP ~ GT+SL+Clutch,data=clean_pitx2_data))
##This makes biological sence, and has the lowest AIC
summary(lm(TVTP ~ GT+SL:Clutch,data=clean_pitx2_data))
##This is technically what is plotted without specifying a fit, https://www.andrew.cmu.edu/user/achoulde/94842/lectures/lecture10/lecture10-94842.html
summary(lm(TVTP ~ GT+SL*Clutch,data=clean_pitx2_data))

##The below uses the fit from the best model
summary(  lm(TVTP ~ GT+SL:Clutch,data=clean_pitx2_data))

tooth_lm<-lm(TVTP ~ GT+SL:Clutch,data=clean_pitx2_data)
tooth_predict <- predict(tooth_lm,interval="confidence")
clean_pitx2_data<-cbind(clean_pitx2_data,predict(tooth_lm,interval="confidence"))
ggplot(clean_pitx2_data, aes(x=SL, y=TVTP,color=GT,linetype=Clutch,shape=Clutch)) + 
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'green','W2'='green','W1'='purple','3737'='red')) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=fit),size=2)

ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_allGT_TVTP_SL.pdf',sep=""))

clean_pitx2_data$tvtp_res <- resid(tooth_lm)
qplot(GT,tvtp_res, geom=c("boxplot","jitter"),color=GT,data = clean_pitx2_data)  +
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'purple','W2'='purple','W1'='purple','3737'='red')) +
  labs(x='GT',y='Corrected Tooth Number',color='Pitx2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_allGT_tvtp_res_box.pdf',sep=""))

sl_clutch_lm <- lm(SL ~ Clutch,data=clean_pitx2_data,na.action = na.exclude)
clean_pitx2_data['sl_res'] <-  resid(sl_clutch_lm,na.action = na.exclude)
qplot(GT,sl_res, geom=c("boxplot","jitter"),color=GT,data = clean_pitx2_data)  +
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'purple','W2'='purple','W1'='purple','3737'='red')) +
  labs(x='GT',y='Corrected_Standard Length',color='Pitx2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_allGT_sl_res_box.pdf',sep=""))

summary(lm(SL~GT+Clutch,data=clean_pitx2_data))
################
##Break it down by GT
############
with(clean_pitx2_data, table(GT, Clutch))
##Let's check out the +37 insertion
clean_pitx2_data_37 <- clean_pitx2_data[clean_pitx2_data$Clutch %in% c('TG1123A','TG1123B','TG1209','TG1264'),]
tooth_lm<-lm(TVTP ~ GT+SL:Clutch,data=clean_pitx2_data_37)
tooth_predict <- predict(tooth_lm,interval="confidence")
clean_pitx2_data_37<-cbind(clean_pitx2_data_37,predict(tooth_lm,interval="confidence"))
ggplot(clean_pitx2_data_37, aes(x=SL, y=TVTP,color=GT,linetype=Clutch,shape=Clutch)) + 
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'purple','W2'='purple','W1'='purple','3737'='red')) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=fit),size=2)
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_37_TVTP_SL.pdf',sep=""))
qplot(GT,SL, geom=c("boxplot","jitter"),color=GT,data = clean_pitx2_data_37)  +
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'purple','W2'='purple','W1'='purple','3737'='red')) +
  labs(x='GT',y='Standard Length',color='Pitx2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 


sl_clutch_lm <- lm(SL ~ Clutch,data=clean_pitx2_data_37,na.action = na.exclude)
clean_pitx2_data_37['sl_res'] <-  resid(sl_clutch_lm,na.action = na.exclude)
qplot(GT,sl_res, geom=c("boxplot","jitter"),color=GT,data = clean_pitx2_data_37)  +
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'purple','W2'='purple','W1'='purple','3737'='red')) +
  labs(x='GT',y='Corrected_Standard Length',color='Pitx2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 

summary(lm(SL~GT+Clutch,data=clean_pitx2_data_37))
TukeyHSD(aov(sl_res~GT,data=clean_pitx2_data_37))

sl_clutch_lm <- lm(TVTP ~ SL:Clutch,data=clean_pitx2_data_37,na.action = na.exclude)
clean_pitx2_data_37['tvtp_res'] <-  resid(sl_clutch_lm,na.action = na.exclude)

qplot(GT,tvtp_res, geom=c("boxplot","jitter"),color=GT,data = clean_pitx2_data_37)  +
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'purple','W2'='purple','W1'='purple','3737'='red')) +
  labs(x='GT',y='Corrected Tooth Number',color='Plod2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 

##Let's check out the -2 Deletion
clean_pitx2_data_2 <- clean_pitx2_data[clean_pitx2_data$Clutch %in% c('TG1168'),]
tooth_lm<-lm(TVTP ~ GT+SL:Clutch,data=clean_pitx2_data_2)
tooth_predict <- predict(tooth_lm,interval="confidence")
clean_pitx2_data_2<-cbind(clean_pitx2_data_2,predict(tooth_lm,interval="confidence"))
ggplot(clean_pitx2_data_2, aes(x=SL, y=TVTP,color=GT,linetype=Clutch,shape=Clutch)) + 
  scale_colour_manual(values = c("WW" = "blue", 'W37' = 'purple','W2'='purple','W1'='purple','3737'='red')) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=fit),size=2)
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_2_TVTP_SL.pdf',sep=""))


##Let's check out the -1 Deletion
clean_pitx2_data_1 <- clean_pitx2_data[clean_pitx2_data$Clutch %in% c('TG1208','TG1227A','TG1227C','TG1322A','TG1322B','TG1323A','TG1323B'),]
tooth_lm<-lm(TVTP ~ GT+SL:Clutch,data=clean_pitx2_data_1)
tooth_predict <- predict(tooth_lm,interval="confidence")
clean_pitx2_data_1<-cbind(clean_pitx2_data_1,predict(tooth_lm,interval="confidence"))
ggplot(clean_pitx2_data_1, aes(x=SL, y=TVTP,color=GT,linetype=Clutch,shape=Clutch)) + 
  scale_colour_manual(values = c("WW" = "royalblue", 'W37' = 'purple','W2'='purple','W1'='magenta','3737'='red')) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=fit),size=2)
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_1_TVTP_SL.pdf',sep=""))

AIC(lm(TVTP ~ GT+SL,data=clean_pitx2_data_1))
AIC(lm(TVTP ~ GT+Clutch,data=clean_pitx2_data_1))
AIC(lm(TVTP ~ GT+SL+Clutch,data=clean_pitx2_data_1))
AIC(lm(TVTP ~ GT+SL:Clutch,data=clean_pitx2_data_1))


sl_lm <- lm(TVTP ~ SL:Clutch,data=clean_pitx2_data_1,na.action = na.exclude)
clean_pitx2_data_1['tvtp_resid'] <-  resid(sl_lm,na.action = na.exclude)
qplot(GT,tvtp_resid, geom=c("boxplot","jitter"),color=GT,data = clean_pitx2_data_1)  +
  scale_colour_manual(values = c("WW" = "royalblue", 'W37' = 'purple','W2'='purple','W1'='magenta','3737'='red')) +
  labs(x='GT',y='Corrected Standard Length',color='Pitx2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_1_box_tvtp_resid.pdf',sep=""))


sl_clutch_lm <- lm(SL ~ Clutch,data=clean_pitx2_data_1,na.action = na.exclude)
clean_pitx2_data_1['sl_res'] <-  resid(sl_clutch_lm,na.action = na.exclude)
qplot(GT,sl_res, geom=c("boxplot","jitter"),color=GT,data = clean_pitx2_data_1)  +
  scale_colour_manual(values = c("WW" = "royalblue", 'W37' = 'purple','W2'='purple','W1'='magenta','3737'='red')) +
  labs(x='GT',y='Corrected Standard Length',color='Pitx2 Genotype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=18)) 
ggsave(paste(Dropbox,'Miller/figures/Mutants/Pitx2/CERC_Pitx2_1_SL.pdf',sep=""))
clean_pitx2_data_1$GT <- as.character(clean_pitx2_data_1$GT)
clean_pitx2_data_1$GT
clean_pitx2_data_1[clean_pitx2_data_1$GT == 'W1',]$GT <- 'WD'
as.factor(clean_pitx2_data_1$GT)
clean_pitx2_data_1$GT <- factor(clean_pitx2_data_1$GT,levels=c('WW','WD','DD'))

ggplot(clean_pitx2_data_1,aes(x=GT,y=sl_res,color=GT)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Length') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Pitx2 \nGenotype",override.aes = list(size=3,linetype=0)))
