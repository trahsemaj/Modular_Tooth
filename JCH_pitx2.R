library(ggplot2)
##PITX2 comb
##NTC Pitx2
#install.packages('ggplot2')
#pitx2_data = read.csv('C:/Users/trahs/Dropbox/Miller/data/Mutants/Pitx2/NTC_PITX2_081817.csv')
pitx2_data = read.csv('/home/james/Dropbox/Miller/data/Mutants/Pitx2/NTC_PITX2_081817.csv')
qplot(SL,TVTP,data=pitx2_data,color=Genotype,shape=Clutch,size=I(3))
pitx2_data <- pitx2_data[!(is.na(pitx2_data$TVTP)),]
pitx2_data <- pitx2_data[!(is.na(pitx2_data$Genotype)),]
backup_options <- options()
options(backup_options)



pitx2_data$Pdiff <- pitx2_data$Pelv_Left - pitx2_data$Pelv_Right

ggplot(pitx2_data, aes(x=SL, y=TVTP,color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(pitx2_data, aes(x=SL, y=Pelv_Left,color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(pitx2_data, aes(x=SL, y=Pelv_Right,color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(pitx2_data, aes(x=SL, y=Pdiff, color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)


AIC(lm(TVTP ~ SL, data=pitx2_data))
AIC(lm(TVTP ~ Clutch, data=pitx2_data))
AIC(lm(TVTP ~ SL + Clutch, data=pitx2_data))
AIC(lm(TVTP ~ SL:Clutch, data=pitx2_data))

summary(lm(TVTP ~ Genotype+SL:Clutch, data=pitx2_data))
summary(aov(lm(TVTP ~ Genotype+SL:Clutch, data=pitx2_data)))
summary(aov(lm(TVTP ~ SL:Clutch+Genotype, data=pitx2_data)))
summary(aov(TVTP ~ SL:Clutch+Genotype, data=pitx2_data))
summary(aov(TVTP ~ Genotype+SL+Clutch, data=pitx2_data))
pitx2_data$TVTP_resid <- resid(lm(TVTP ~ SL:Clutch, data=pitx2_data,na.action=na.exclude))
summary(aov(TVTP_resid ~ Genotype, data=pitx2_data))
TukeyHSD(aov(TVTP_resid ~ Genotype, data=pitx2_data))

options(contrasts = c('contr.sum','contr.poly'))
pitx2_lm <- lm(TVTP ~ SL:Clutch + Genotype, data=pitx2_data)
drop1(pitx2_lm, .~., test='F')
options(contrasts = c('contr.sum','contr.poly'))
pitx2_data$Genotype


options(contrasts = c('contr.sum','contr.poly'))


pitx2_data$Genotype <- factor(pitx2_data$Genotype,levels=c('WW','W1','W2','W37','3737'))
summary(lm(SL~Genotype+Clutch,data=pitx2_data))
pitx2_sl_lm <- lm(SL ~ Clutch + Genotype, data=pitx2_data)
drop1(pitx2_sl_lm, .~., test='F')
summary(pitx2_sl_lm)
options(backup_options)

plot(TVTP_resid ~ Genotype, data=pitx2_data)
points(TVTP_resid ~ Genotype, data=pitx2_data)

pitx2_data$Simple_GT <- NA
pitx2_data[which(pitx2_data$Genotype == 'WW'),]$Simple_GT <- 'WT'
pitx2_data[which(pitx2_data$Genotype == 'W1'),]$Simple_GT <- 'Mut'
pitx2_data[which(pitx2_data$Genotype == 'W2'),]$Simple_GT <- 'Mut'
pitx2_data[which(pitx2_data$Genotype == 'W37'),]$Simple_GT <- 'Mut'
pitx2_data[which(pitx2_data$Genotype == '3737'),]$Simple_GT <- 'Mut'
pitx2_data$Simple_GT
pitx2_data$Simple_GT <- factor(pitx2_data$Simple_GT,levels=c('WT','Mut'))
pitx2_data$Simple_GT


options(contrasts = c('contr.sum','contr.poly'))
summary(lm(TVTP ~ SL:Clutch + Simple_GT,data=pitx2_data))
test_lm <- lm(TVTP ~ SL:Clutch + Simple_GT, data=pitx2_data)
drop1(test_lm, .~., test='F')

ggplot(pitx2_data, aes(x=SL, y=TVTP,color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(pitx2_data, aes(x=SL, y=TVTP,color=Genotype)) + 
  geom_point(size=3) + 
  theme_linedraw() +
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple','W1' = 'purple','W2' = 'orange','W37' = 'green','3737' = 'red')) +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=14)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(pitx2_data, aes(x=Genotype, y=TVTP_resid,color=Genotype)) + 
  geom_boxplot(outlier.size=0) +
  geom_jitter() +
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple','W1' = 'purple','W2' = 'orange','W37' = 'green','3737' = 'red')) +
  labs(x='Clutch',y='Corrected Tooth Number',color='Pitx2 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14))

ggplot(pitx2_data, aes(x=Clutch, y=SL,color=Genotype)) + 
  geom_boxplot(outlier.size=0) +
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple','W1' = 'purple','W2' = 'orange','W37' = 'green','3737' = 'red')) +
  labs(x='Clutch',y='SL',color='Pitx2 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 


ggplot(pitx2_data, aes(x=Clutch, y=TVTP_resid,color=Genotype)) + 
  geom_boxplot(outlier.size=0) +
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple','W1' = 'purple','W2' = 'orange','W37' = 'green','3737' = 'red')) +
  labs(x='Clutch',y='Corrected TVTP',color='Pitx2 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 


##TG1227 includes SAM data

TG1227 <- pitx2_data[which(pitx2_data$Clutch %in% c('TG1227A','TG1227C','TG1208','TG1322A','TG1322B','TG1323A','TG1323B')),]
#TG1227 <- pitx2_data[which(pitx2_data$Clutch %in% c('TG1227A','TG1227C','TG1208')),]
AIC(lm(TVTP ~ SL, data=TG1227))
AIC(lm(TVTP ~ Clutch, data=TG1227))
AIC(lm(TVTP ~ SL + Clutch, data=TG1227))
AIC(lm(TVTP ~ SL:Clutch, data=TG1227))

AIC(lm(TVTP ~ SL:Clutch+Genotype, data=TG1227))
AIC(lm(TVTP ~ Genotype*Clutch*SL, data=TG1227))


summary(lm(TVTP ~ SL + Clutch, data=TG1227))
summary(lm(TVTP ~ SL:Clutch, data=TG1227))
summary(aov(lm(TVTP ~ Genotype+Clutch+SL, data=TG1227)))
summary(aov(lm(TVTP ~ Clutch+SL+Genotype, data=TG1227)))
summary(aov(lm(TVTP ~ Clutch:SL+Genotype, data=TG1227)))
summary(aov(lm(TVTP ~ Genotype+Clutch:SL, data=TG1227)))
summary(aov(lm(TVTP ~ Genotype, data=TG1227)))

TG1227$TVTP_resid <- resid(lm(TVTP ~ SL:Clutch, data=TG1227,na.action=na.exclude))
wilcox.test(TG1227[which(TG1227$Genotype == 'WW'),]$TVTP_resid,
TG1227[which(TG1227$Genotype == 'W1'),]$TVTP_resid)


summary(aov(lm(TVTP_resid ~ Genotype, data=TG1227)))
t.test(TG1227[which(TG1227$Genotype == 'WW'),]$TVTP_resid,
            TG1227[which(TG1227$Genotype == 'W1'),]$TVTP_resid)

summary(aov(lm(TVTP ~ SL:Clutch + Genotype, data=TG1227)))
summary(aov(lm(TVTP ~ Genotype+SL:Clutch, data=TG1227)))
summary(aov(lm(TVTP ~ SL:Clutch+Genotype, data=TG1227)))
summary(aov(TVTP_resid ~ Genotype, data=TG1227))
summary(aov(TVTP ~ Genotype + Clutch+SL, data=TG1227))
summary(aov(TVTP ~ SL+Clutch+Genotype , data=TG1227))


AIC(aov(TVTP ~ Genotype+Clutch+SL, data=TG1227))
AIC(aov(TVTP ~ SL:Clutch+Genotype, data=TG1227))
aggregate(TVTP~Genotype+Clutch,data=TG1227,length)$TVTP

##Type III anova on the best LM
options(contrasts = c('contr.sum','contr.poly'))
TG1227_lm <- lm(TVTP ~ SL:Clutch + Genotype, data=TG1227)
drop1(TG1227_lm, .~., test='F')
options(backup_options)

options(contrasts = c('contr.sum','contr.poly'))
TG1227_sl_lm <- lm(SL ~ Clutch + Genotype, data=TG1227)
drop1(TG1227_sl_lm, .~., test='F')
options(backup_options)



TG1227$Genotype <- factor(TG1227$Genotype,levels=c('WW','W1'))
table(TG1227$Genotype)
ggplot(TG1227, aes(x=Clutch, y=TVTP,color=Genotype)) + 
  geom_boxplot(outlier.size=0) +
  geom_jitter(position=position_dodge(0.8)) +
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple')) +
  labs(x='Clutch',y='Corrected Tooth Number',color='Pitx2 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 

summary(lm(TVTP ~ SL:Clutch + Genotype, data=TG1227))
summary(lm(TVTP ~ Genotype + Clutch, data=TG1227))
summary(lm(TVTP ~ Genotype + Clutch+SL, data=TG1227))
summary(lm(TVTP ~ SL:Clutch + Genotype, data=TG1227))



summary(glm(TVTP ~ SL:Clutch + Genotype, data=TG1227))
summary(aov(glm(TVTP ~ SL+Clutch + Genotype, data=TG1227)))

ggplot(TG1227, aes(x=SL, y=TVTP,color=Genotype)) + 
  geom_point(size=2) + 
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple')) +
  labs(x='SL',y='Tooth Number',color='Pitx2 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) +
  geom_smooth(method=lm,se=FALSE)


#ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/Pitx2_TG1208_1227_TVTP.pdf',width=6,height=4.5,units='in')


ggplot(TG1227, aes(x=Clutch, y=SL,color=Genotype)) + 
  geom_boxplot(outlier.size=NA) +
  geom_jitter(position=position_dodge(0.8)) +
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple')) +
  labs(x='Clutch',y='Standard Length',color='Pitx2 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/Pitx2_TG1208_1227_SL.pdf',width=6,height=4.5,units='in')


ggplot(TG1227, aes(x=SL, y=TVTP,color=Genotype)) + 
  scale_colour_manual(values = c("WW" = "blue", 'W1' = 'purple')) +
  geom_point(size=1) + 
  theme_linedraw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/Pitx2_TG1208_1227_TVTP_vs_SL.pdf',width=6,height=4.5,units='in')


ggplot(TG1227, aes(x=SL, y=Pelv_Left,color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(TG1227, aes(x=SL, y=Pelv_Right,color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(TG1227, aes(x=SL, y=Pdiff,color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)



ggplot(pitx2_data, aes(x=SL, y=TVTP,color=Genotype)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)


AIC(lm(TVTP ~ SL, data=TG1227))
AIC(lm(TVTP ~ Clutch, data=TG1227))
AIC(lm(TVTP ~ SL + Clutch, data=TG1227))
AIC(lm(TVTP ~ SL*Clutch, data=TG1227))
AIC(lm(TVTP ~ SL:Clutch, data=TG1227))


AIC(lm(Pelv_Right ~ SL, data=TG1227))


TG1227$TVTP_resid <- resid(lm(TVTP ~ SL *Clutch, data=TG1227,na.action=na.exclude))
TG1227$PL_resid <- resid(lm(Pelv_Left ~ SL , data=TG1227,na.action=na.exclude))
TG1227$PR_resid <- resid(lm(Pelv_Right ~ SL , data=TG1227,na.action=na.exclude))



plot(TVTP_resid ~ Genotype, data=TG1227)
points(TVTP_resid ~ Genotype, data=TG1227)
summary(aov(TVTP_resid ~ Genotype, data=TG1227))
wilcox.test(TG1227[which(TG1227$Genotype == 'WW'),]$TVTP_resid,
TG1227[which(TG1227$Genotype == 'W1'),]$TVTP_resid)


plot(PL_resid ~ Genotype, data=TG1227)
points(PL_resid ~ Genotype, data=TG1227)
summary(aov(PL_resid ~ Genotype, data=TG1227))


plot(PR_resid ~ Genotype, data=TG1227)
points(PR_resid ~ Genotype, data=TG1227)
summary(aov(PR_resid ~ Genotype, data=TG1227))

summary(aov(Pdiff ~ Genotype, data=TG1227))
##SAM Pitx2 stuff below here
tooth_data_37 = read.csv('/home/james/Dropbox/Miller/data/Mutants/Pitx2/Siegen_Pitx2_TG1209_TG1123.csv')
tooth_data_1 = read.csv('/home/james/Dropbox/Miller/data/Mutants/Pitx2/Siegen_Pitx2_TG1208.csv')

tooth_data_37$Class <- factor(tooth_data_37$Class, levels = levels(tooth_data_37$Class)[c(2,1,3)])
table(tooth_data_1$Tank )
table(tooth_data_37[which(tooth_data_37$Tank == 'C'),]$Class)
plot(TVTP ~ SL,data=tooth_data_37)
plot(TVTP ~ SL,data=tooth_data_1)

qplot(SL,TVTP,data=tooth_data_37,color=Class,shape=Tank,size=I(3))
qplot(SL,TVTP,data=tooth_data_1,color=Class,size=I(3))

ggplot(tooth_data_37, aes(x=SL, y=TVTP,color=Class,shape=Tank)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

ggplot(tooth_data_1, aes(x=SL, y=TVTP,color=Class)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

summary(lm(TVTP ~ SL,data=tooth_data_37))
summary(lm(TVTP ~ Tank,data=tooth_data_37))
summary(lm(TVTP ~ SL+Tank,data=tooth_data_37))
summary(lm(TVTP ~ SL*Tank,data=tooth_data_37))

summary(lm(TVTP ~ SL+Tank+Class,data=tooth_data_37,na.action=na.exclude))
lm_37_sl_tank <- lm(TVTP ~ SL+Tank,data=tooth_data_37,na.action=na.exclude)
tooth_data_37['tvtp_res'] <- resid(lm_37_sl_tank,na.action=na.exclude)
plot(tvtp_res ~ Class,data=tooth_data_37)
summary(aov(tvtp_res ~ Class,data=tooth_data_37))


qplot(Class,tvtp_res, geom=c("boxplot","jitter"),color=Class,data = tooth_data_37)  + theme(text = element_text(size=20))
barplot(table(tooth_data_37[which(tooth_data_37$Tank == 'C'),]$Class))
summary(lm(TVTP ~ SL,data=tooth_data_1))
toplot <- tooth_data_37[which(tooth_data_37$Tank == 'C'),]
qplot(Class,SL, color=Class,geom=c("boxplot","jitter"),data=toplot)  + theme(text = element_text(size=20))


summary(lm(TVTP ~ SL+Class,data=tooth_data_1,na.action=na.exclude))
lm_1_sl <- lm(TVTP ~ SL,data=tooth_data_1,na.action=na.exclude)
tooth_data_1['tvtp_res'] <- resid(lm_1_sl ,na.action=na.exclude)
plot(tvtp_res ~ Class,data=tooth_data_1)
qplot(Class,tvtp_res, geom=c("boxplot","jitter"),color=Class,data = tooth_data_1) + theme(text = element_text(size=20))
summary(aov(tvtp_res ~ Class,data=tooth_data_1))
summary(aov(SL ~ Class,data=tooth_data_1))
qplot(Class,SL, geom=c("boxplot","jitter"),color=Class,data = tooth_data_1)  + theme(text = element_text(size=20))
