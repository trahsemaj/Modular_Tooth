library(ggplot2)



plod2_data <- read.csv('/home/james/Dropbox/Miller/data/Mutants/Plod2/Plod2_TG1391_oral.csv')
plod2_data[is.na(plod2_data$Sex),]
plod2_data <- plod2_data[!is.na(plod2_data$GT),]
plod2_data <- plod2_data[!is.na(plod2_data$SL),]
plod2_data <- plod2_data[!is.na(plod2_data$TVTP),]
plod2_data <- plod2_data[!is.na(plod2_data$Sex),]
plod2_data$GT_2 <- 'wt'
plod2_data[which(plod2_data$GT=='DD'),]$GT_2 <- 'mut'
plod2_data$GT_2 <- as.factor(plod2_data$GT_2)
plod2_data$GT_2 


plod2_data_wt <- plod2_data[which(plod2_data$GT == 'WW'),]
ggplot(plod2_data_wt, aes(x=SL, y=TotalDentary,color=Sex)) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)
ggplot(plod2_data_wt, aes(x=SL, y=TotalPreMax,color=Sex)) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(plod2_data, aes(x=SL, y=TotalDentary,color=GT,shape=Sex,linetype=Sex)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(plod2_data, aes(x=SL, y=TotalPreMax,color=GT,shape=Sex,linetype=Sex)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(plod2_data, aes(x=SL, y=TDTP1,color=GT,shape=Sex,linetype=Sex)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(plod2_data, aes(x=SL, y=TDTP2,color=GT,shape=Sex,linetype=Sex)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

summary(lm(TotalDentary~ Sex + SL + GT,data=plod2_data))
summary(lm(TotalPreMax~ Sex + SL + GT,data=plod2_data))

summary(lm(TotalPreMax~ Sex + SL + GT_2,data=plod2_data))

summary(lm(TVTP~ Sex + SL + GT,data=plod2_data))
summary(lm(TVTP~ Sex + SL + GT_2,data=plod2_data))

summary(lm(TDTP1~ Sex + SL + GT,data=plod2_data))
summary(lm(TDTP1~ Sex + SL + GT_2,data=plod2_data))

summary(lm(TDTP2~ Sex + SL + GT,data=plod2_data))
summary(lm(TDTP2~ Sex + SL + GT_2,data=plod2_data))

plod2_data[which(plod2_data$GT=='DD' & plod2_data$TotalPreMax > 40),]

ggplot(plod2_data, aes(x=Dentary, y=Pre.Max,color=GT)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))

cor(plod2_data_wt$Dentary,plod2_data_wt$Pre.Max,method='spearman')
cor(plod2_data_wt$Dentary,plod2_data_wt$Pre.Max,method='spearman')


ggplot(plod2_data, aes(x=Dentary, y=TVTP,color=GT)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))

ggplot(plod2_data, aes(x=Pre.Max, y=TVTP,color=GT)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))


cor(plod2_data$TVTP,plod2_data$Pre.Max,method='spearman')
cor(plod2_data_wt$TVTP,plod2_data_wt$Pre.Max,method='spearman')


cor(plod2_data$TVTP,plod2_data$Dentary,method='spearman')
cor(plod2_data_wt$TVTP,plod2_data_wt$Dentary,method='spearman')

ggplot(plod2_data, aes(x=SL, y=TVTP,color=GT)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))

ggplot(plod2_data, aes(x=SL, y=Dentary,color=GT)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))


ggplot(plod2_data, aes(x=SL, y=Pre.Max,color=GT)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))


##fish 26
plod2_data[which(plod2_data$GT == 'DD' & plod2_data$Pre.Max > 22),]

summary(lm(TVTP ~ GT + SL, data=plod2_data))
summary(lm(Dentary ~ GT + SL, data=plod2_data))
summary(lm(Pre.Max ~ GT + SL, data=plod2_data))

