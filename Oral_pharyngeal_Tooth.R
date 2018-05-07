library(ggplot2)
library(reshape2)
backup_options <- options()

Folder = '/home/james/Dropbox/Miller/data/Tooth_Counts/'
Folder = 'C:/Users/trahs/Dropbox/Miller/data/Tooth_Counts/'


total_tooth_data <- read.csv(paste(Folder,'Filtered_Pop_Tooth_Counts_oral_filtered.csv',sep=''))
sixpop_tooth_data <- total_tooth_data[which(total_tooth_data$Pop %in% c('JAMA','LITC','RABS','FTC','CERC','PAXB')),]
sixpop_tooth_data$Pop <- factor(sixpop_tooth_data$Pop,levels=c('RABS','JAMA','LITC','FTC','CERC','PAXB'))
sixpop_tooth_data <- sixpop_tooth_data[!is.na(sixpop_tooth_data$TVTP),]
sixpop_tooth_data <- sixpop_tooth_data[!is.na(sixpop_tooth_data$TDTP1),]
sixpop_tooth_data <- sixpop_tooth_data[!is.na(sixpop_tooth_data$TDTP2),]
sixpop_tooth_data <- sixpop_tooth_data[sixpop_tooth_data$SL>30,]
sixpop_tooth_data$Pop_Type <- 'Marine'
sixpop_tooth_data$Pop_Type
sixpop_tooth_data[which(sixpop_tooth_data$Pop=='CERC'),]$Pop_Type = 'CERC'
sixpop_tooth_data[which(sixpop_tooth_data$Pop=='PAXB'),]$Pop_Type = 'PAXB'
sixpop_tooth_data[which(sixpop_tooth_data$Pop=='FTC'),]$Pop_Type = 'FTC'
sixpop_tooth_data$Pop_Type <- factor(sixpop_tooth_data$Pop_Type,levels=c('Marine','FTC','CERC','PAXB'))
summary(lm(TVTP~SL+Pop,data=sixpop_tooth_data))
summary(lm(TDTP1~SL +Pop+Sex,data=sixpop_tooth_data))
summary(lm(TDTP2~SL +Pop,data=sixpop_tooth_data))

summary(lm(TVTP~SL+Pop_Type,data=sixpop_tooth_data))
summary(lm(TDTP1~SL +Pop_Type,data=sixpop_tooth_data))
summary(lm(TDTP2~SL +Pop_Type,data=sixpop_tooth_data))



sixpop_tooth_data$TVTP_resid <- resid(lm(TVTP~SL,data=sixpop_tooth_data))
sixpop_tooth_data$TDTP1_resid <- resid(lm(TDTP1~SL,data=sixpop_tooth_data))
sixpop_tooth_data$TDTP2_resid <- resid(lm(TDTP2~SL,data=sixpop_tooth_data))

Outfolder = '/home/james/Dropbox/Miller/figures/Multipop_intron4/Tooth_Counts/'
Outfolder = 'C:/Users/trahs/Dropbox/Miller/figures/Multipop_intron4/Tooth_Counts/'

tvtp_lm<-lm(TVTP ~ SL+Pop,data=sixpop_tooth_data)
tvtp_predict <- predict(tvtp_lm,interval="confidence")
colnames(tvtp_predict) <- c('tvtp_fit','tvtp_lwr','tvtp_upr')
sixpop_tooth_data<-cbind(sixpop_tooth_data,tvtp_predict)
ggplot(sixpop_tooth_data, aes(x=SL, y=TVTP,color=Pop)) + 
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='steelblue','JAMA'='darkred')) +
  geom_point(size=1.5) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(x=SL,y=tvtp_fit,color=Pop),size=1,inherit.aes=FALSE)
ggsave(paste(Outfolder,'6pop_SL_TVTP.pdf'))

tdtp1_lm<-lm(TDTP1 ~ SL+Pop,data=sixpop_tooth_data)
tdtp1_predict <- predict(tdtp1_lm,interval="confidence")
colnames(tdtp1_predict) <- c('tdtp1_fit','tdtp1_lwr','tdtp1_upr')
sixpop_tooth_data<-cbind(sixpop_tooth_data,tdtp1_predict)
ggplot(sixpop_tooth_data[!is.na(sixpop_tooth_data$TDTP1),], aes(x=SL, y=TDTP1,color=Pop)) + 
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='steelblue','JAMA'='darkred')) +
  geom_point(size=1.5) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(x=SL,y=tdtp1_fit,color=Pop),size=1,inherit.aes=FALSE)
ggsave(paste(Outfolder,'6pop_SL_TDTP1.pdf'))

tdtp2_lm<-lm(TDTP2 ~ SL+Pop,data=sixpop_tooth_data)
tdtp2_predict <- predict(tdtp2_lm,interval="confidence")
colnames(tdtp2_predict) <- c('tdtp2_fit','tdtp2_lwr','tdtp2_upr')
sixpop_tooth_data<-cbind(sixpop_tooth_data,tdtp2_predict)
ggplot(sixpop_tooth_data, aes(x=SL, y=TDTP2,color=Pop)) + 
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='steelblue','JAMA'='darkred')) +
  geom_point(size=1.5) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(x=SL,y=tdtp2_fit,color=Pop),size=1,inherit.aes=FALSE)
ggsave(paste(Outfolder,'6pop_SL_TDTP2.pdf'))

ggplot(sixpop_tooth_data,aes(x=Pop,y=TVTP_resid,color=Pop)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=.5) +
  geom_jitter(size=.5) +
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='steelblue','JAMA'='darkred')) +
  theme_bw() +
  labs(y='Corrected VTP') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30))
ggsave(paste(Outfolder,'6pop_TVTP_res_box.pdf'))


ggplot(sixpop_tooth_data,aes(x=Pop,y=TDTP1_resid,color=Pop)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=.5) +
  geom_jitter(size=.5) +
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='steelblue','JAMA'='darkred')) +
  theme_bw() +
  labs(y='Corrected DTP1') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30))
ggsave(paste(Outfolder,'6pop_TDTP1_res_box.pdf'))


ggplot(sixpop_tooth_data,aes(x=Pop,y=TDTP2_resid,color=Pop)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=.5) +
  geom_jitter(size=.5) +
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='steelblue','JAMA'='darkred')) +
  theme_bw() +
  labs(y='Corrected DTP2') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30))
ggsave(paste(Outfolder,'6pop_TDTP2_res_box.pdf'))

fourpop_tooth_data <- total_tooth_data[which(total_tooth_data$Pop %in% c('LITC','RABS','CERC','PAXB')),]
fourpop_tooth_data <- fourpop_tooth_data[!is.na(fourpop_tooth_data$TVTP),]
fourpop_tooth_data <- fourpop_tooth_data[!is.na(fourpop_tooth_data$TDTP1),]
fourpop_tooth_data <- fourpop_tooth_data[!is.na(fourpop_tooth_data$TDTP2),]
fourpop_tooth_data <- fourpop_tooth_data[!is.na(fourpop_tooth_data$Premax),]
fourpop_tooth_data <- fourpop_tooth_data[!is.na(fourpop_tooth_data$Dentary),]
fourpop_tooth_data <- fourpop_tooth_data[fourpop_tooth_data$SL>30,]



summary(lm(TVTP~SL + Sex + Pop,data=fourpop_tooth_data))
summary(lm(TDTP1~SL + Sex + Pop,data=fourpop_tooth_data))
summary(lm(TDTP2~SL + Sex + Pop,data=fourpop_tooth_data))
summary(lm(Premax~SL + Sex + Pop,data=fourpop_tooth_data))
summary(lm(Dentary~SL + Sex + Pop,data=fourpop_tooth_data))

##Just males
summary(lm(TVTP~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'M'),]))
summary(lm(TDTP1~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'M'),]))
summary(lm(TDTP2~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'M'),]))
summary(lm(Premax~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'M'),]))
summary(lm(Dentary~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'M'),]))

##Just females
summary(lm(TVTP~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'F'),]))
summary(lm(TDTP1~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'F'),]))
summary(lm(TDTP2~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'F'),]))
summary(lm(Premax~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'F'),]))
summary(lm(Dentary~SL + Pop,data=fourpop_tooth_data[which(fourpop_tooth_data$Sex == 'F'),]))

oral_tooth_data <- fourpop_tooth_data

oral_tooth_data$TOT <- oral_tooth_data$Premax + oral_tooth_data$Dentary

oral_tooth_data$Pop <- factor(oral_tooth_data$Pop,levels=c('RABS','LITC','CERC','PAXB'))
oral_tooth_data$Pop_Type <- 'Marine'
oral_tooth_data$Pop_Type
oral_tooth_data[which(oral_tooth_data$Pop=='CERC'),]$Pop_Type = 'CERC'
oral_tooth_data[which(oral_tooth_data$Pop=='PAXB'),]$Pop_Type = 'PAXB'
oral_tooth_data$Pop_Type <- factor(oral_tooth_data$Pop_Type,levels=c('Marine','CERC','PAXB'))
summary(lm(TVTP~SL*Sex,data=oral_tooth_data))
summary(lm(TDTP1~SL*Sex,data=oral_tooth_data))
summary(lm(TDTP2~SL * Sex,data=oral_tooth_data))
summary(lm(Premax~SL +Sex*Pop,data=oral_tooth_data))
summary(lm(Dentary~SL+Sex*Pop,data=oral_tooth_data))


oral_tooth_data$Premax_resid <- resid(lm(Premax~SL+Sex,data=oral_tooth_data))
oral_tooth_data$Dentary_resid <- resid(lm(Dentary~SL+Sex,data=oral_tooth_data))
oral_tooth_data$TOT_resid <- resid(lm(TOT~SL+Sex,data=oral_tooth_data))

summary(lm(Premax~SL + Sex+Pop_Type,data=oral_tooth_data))
summary(lm(Dentary~SL + Sex+Pop_Type,data=oral_tooth_data))
summary(lm(TOT~SL + Sex+Pop_Type,data=oral_tooth_data))

summary(lm(Premax~SL + Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'M'),]))
summary(lm(Dentary~SL + Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'M'),]))
summary(lm(TOT~SL + Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'M'),]))


summary(lm(Premax~SL + Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'F'),]))
summary(lm(Dentary~SL + Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'F'),]))
summary(lm(TOT~SL + Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'F'),]))

premax_aov <- aov(Premax_resid ~ Pop,data=oral_tooth_data)
summary(premax_aov)
TukeyHSD(premax_aov)

dentary_aov <- aov(Dentary_resid ~ Pop,data=oral_tooth_data)
summary(dentary_aov)
TukeyHSD(dentary_aov)

tot_aov <- aov(TOT_resid ~ Pop,data=oral_tooth_data)
summary(tot_aov)
TukeyHSD(tot_aov)

premax_aov <- aov(Premax_resid ~ Pop_Type,data=oral_tooth_data)
summary(premax_aov)
TukeyHSD(premax_aov)

premax_aov <- aov(Premax_resid ~ Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'M'),])
summary(premax_aov)
TukeyHSD(premax_aov)

premax_aov <- aov(Premax_resid ~ Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'F'),])
summary(premax_aov)
TukeyHSD(premax_aov)

dentary_aov <- aov(Dentary_resid ~ Pop_Type,data=oral_tooth_data)
summary(dentary_aov)
TukeyHSD(dentary_aov)

dentary_aov <- aov(Dentary_resid ~ Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'M'),])
summary(dentary_aov)
TukeyHSD(dentary_aov)

dentary_aov <- aov(Dentary_resid ~ Pop_Type,data=oral_tooth_data[which(oral_tooth_data$Sex == 'F'),])
summary(dentary_aov)
TukeyHSD(dentary_aov)

tot_aov <- aov(TOT_resid ~ Pop_Type,data=oral_tooth_data)
summary(tot_aov)
TukeyHSD(tot_aov)


Outfolder = '/home/james/Dropbox/Miller/figures/Multipop_intron4/Tooth_Counts/'
Outfolder = 'C:/Users/trahs/Dropbox/Miller/figures/Multipop_intron4/Tooth_Counts/'
ggplot(oral_tooth_data,aes(x=Pop,y=TOT_resid,color=Pop)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=.5) +
  geom_jitter(size=.5) +
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='skyblue','JAMA'='purple')) +
  theme_bw() +
  labs(y='Corrected TOT') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30))
ggsave(paste(Outfolder,'4pop_TOT_res_box.pdf'))


ggplot(oral_tooth_data,aes(x=Pop,y=Premax_resid,color=Pop)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=.5) +
  geom_jitter(size=.5) +
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='skyblue','JAMA'='purple')) +
  theme_bw() +
  labs(y='Corrected Premax') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30))
ggsave(paste(Outfolder,'4pop_Premax_res_box.pdf'))


ggplot(oral_tooth_data,aes(x=Pop,y=Dentary_resid,color=Pop)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=.5) +
  geom_jitter(size=.5) +
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red','FTC'='skyblue','JAMA'='purple')) +
  theme_bw() +
  labs(y='Corrected Dentary') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30))
ggsave(paste(Outfolder,'4pop_Dentary_res_box.pdf'))


den_lm<-lm(Dentary ~ Sex+SL,data=oral_tooth_data)
den_predict <- predict(den_lm,interval="confidence")
colnames(den_predict) <- c('den_fit','den_lwr','den_upr')
oral_tooth_data<-cbind(oral_tooth_data,den_predict)
ggplot(oral_tooth_data, aes(x=SL, y=Dentary,color=Pop,shape=Sex,linetype=Sex)) + 
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red')) +
  geom_point(size=1.5) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(x=SL,y=den_fit,linetype=Sex),size=1,inherit.aes=FALSE)
ggsave(paste(Outfolder,'4pop_SL_Dentary.pdf'))

pm_lm<-lm(Premax ~ Sex+SL,data=oral_tooth_data)
pm_predict <- predict(pm_lm,interval="confidence")
colnames(pm_predict) <- c('pm_fit','pm_lwr','pm_upr')
oral_tooth_data<-cbind(oral_tooth_data,pm_predict)
ggplot(oral_tooth_data, aes(x=SL, y=Premax,color=Pop,shape=Sex,linetype=Sex)) + 
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red')) +
  geom_point(size=1.5) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(x=SL,y=pm_fit,linetype=Sex),size=1,inherit.aes=FALSE)
ggsave(paste(Outfolder,'4pop_SL_Premax.pdf'))

tot_lm<-lm(TOT ~ Sex+SL,data=oral_tooth_data)
tot_predict <- predict(tot_lm,interval="confidence")
colnames(tot_predict) <- c('tot_fit','tot_lwr','tot_upr')
oral_tooth_data<-cbind(oral_tooth_data,tot_predict)
ggplot(oral_tooth_data, aes(x=SL, y=TOT,color=Pop,shape=Sex,linetype=Sex)) + 
  scale_colour_manual(values = c("CERC" = "forestgreen", 'PAXB' = 'blue', "RABS" = "purple",'LITC'='red')) +
  geom_point(size=1.5) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(x=SL,y=tot_fit,linetype=Sex),size=1,inherit.aes=FALSE)
ggsave(paste(Outfolder,'4pop_SL_TOT.pdf'))


options(contrasts = c('contr.sum','contr.poly'))
tvtp_lm <- lm(TVTP~SL+Sex+Pop,data=oral_tooth_data)
drop1(tvtp_lm, .~., test='F')
summary(tvtp_lm)
options(backup_options)

options(contrasts = c('contr.sum','contr.poly'))
tdtp1_lm <- lm(TDTP1~SL+Sex+Pop,data=oral_tooth_data)
drop1(tdtp1_lm, .~., test='F')
summary(tdtp1_lm)
options(backup_options)


options(contrasts = c('contr.sum','contr.poly'))
tdtp2_lm <- lm(TDTP2~SL+Sex+Pop,data=oral_tooth_data)
drop1(tdtp2_lm, .~., test='F')
summary(tvtp_lm)
options(backup_options)

options(contrasts = c('contr.sum','contr.poly'))
den_lm <- lm(Dentary~SL+Sex+Pop,data=oral_tooth_data)
drop1(den_lm, .~., test='F')
summary(den_lm)
options(backup_options)

options(contrasts = c('contr.sum','contr.poly'))
pm_lm <- lm(Premax~SL+Sex+Pop,data=oral_tooth_data)
drop1(pm_lm, .~., test='F')
summary(pm_lm)
options(backup_options)

summary(lm(TVTP~SL+Sex+Pop,data=oral_tooth_data))
summary(lm(TDTP1~SL+Sex+Pop,data=oral_tooth_data))
summary(lm(TDTP2~SL+Sex+Pop,data=oral_tooth_data))
summary(lm(Premax~SL+Sex+Pop,data=oral_tooth_data))
summary(lm(Dentary~SL+Sex+Pop,data=oral_tooth_data))
summary(lm(TOT~SL+Sex+Pop,data=oral_tooth_data))




oral_tooth_data$TVTP_resid <- resid(lm(TVTP~SL+Sex,data=oral_tooth_data))
oral_tooth_data$TDTP1_resid <- resid(lm(TDTP1~SL+Sex,data=oral_tooth_data))
oral_tooth_data$TDTP2_resid <- resid(lm(TDTP2~SL+Sex,data=oral_tooth_data))
oral_tooth_data$Dentary_resid <- resid(lm(Premax~SL+Sex,data=oral_tooth_data))
oral_tooth_data$Premax_resid <- resid(lm(Dentary~SL+Sex,data=oral_tooth_data))

tooth_cols = c('TVTP_resid','TDTP1_resid','TDTP2_resid','Dentary_resid','Premax_resid')
tooth_cols = c('TVTP','TDTP1','TDTP2','Dentary','Premax')
cor.test(oral_tooth_data$TVTP_resid,oral_tooth_data$Dentary_resid)
oral_tooth_cor <- cor(oral_tooth_data[tooth_cols])
oral_tooth_cor <- round(oral_tooth_cor,2)
dist <- as.dist(1-oral_tooth_cor)
hclust_order <- hclust(dist)
hclust_order 
oral_tooth_cor <- oral_tooth_cor[hclust_order$order,hclust_order$order]
oral_tooth_cor[lower.tri(oral_tooth_cor)] <- NA
oral_tooth_cor_melt <- melt(oral_tooth_cor, na.rm = TRUE)
oral_tooth_cor_melt
ggheatmap <- ggplot(oral_tooth_cor_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = .5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
#Outfolder = '/home/james/Dropbox/Miller/figures/Multipop_intron4/Tooth_Counts/'
Outfolder = 'C:/Users/trahs/Dropbox/Miller/figures/Multipop_intron4/Tooth_Counts/'
ggsave(paste(Outfolder,'4pop_tooth_cor_heatmap_slcorrect.pdf'))



F2_data_file = 'C:/Users/trahs/Dropbox/Miller/data/QTL_Data/CxL_Mapping/SAM_CxL_phenos_raw.csv'
F2_data_file = '/home/james/Dropbox/Miller/data/QTL_Data/CxL_Mapping/SAM_CxL_phenos_raw.csv'
cxl_F2_data <- read.csv(F2_data_file)
cxl_F2_data <- cxl_F2_data[!(is.na(cxl_F2_data$TVTP)),]
cxl_F2_data <- cxl_F2_data[!(is.na(cxl_F2_data$TDTP1)),]
cxl_F2_data <- cxl_F2_data[!(is.na(cxl_F2_data$TDTP2)),]
cxl_F2_data <- cxl_F2_data[!(is.na(cxl_F2_data$Premax)),]
cxl_F2_data <- cxl_F2_data[!(is.na(cxl_F2_data$Dentary)),]

ggplot(cxl_F2_data, aes(x=SL, y=TVTP,shape=Sex,linetype=Sex)) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE) 

ggplot(cxl_F2_data, aes(x=SL, y=TDTP1,shape=Sex,linetype=Sex)) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE) 

ggplot(cxl_F2_data, aes(x=SL, y=TDTP2,shape=Sex,linetype=Sex)) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE) 

ggplot(cxl_F2_data, aes(x=SL, y=Premax,shape=Sex,linetype=Sex)) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE) 

ggplot(cxl_F2_data, aes(x=SL, y=Dentary,shape=Sex,linetype=Sex)) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE) 

ggplot(cxl_F2_data, aes(x=SL, y=Total.Oral.Teeth,shape=Sex,linetype=Sex)) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE) 

summary(lm(TVTP~SL+Sex,data=cxl_F2_data))
summary(lm(TDTP1~SL+Sex,data=cxl_F2_data))
summary(lm(TDTP2~SL+Sex,data=cxl_F2_data))
summary(lm(Premax~SL+Sex,data=cxl_F2_data))
summary(lm(Dentary~SL+Sex,data=cxl_F2_data))

cxl_F2_data$TVTP_resid <- resid(lm(TVTP~SL+Sex,data=cxl_F2_data))
cxl_F2_data$TDTP1_resid <- resid(lm(TDTP1~SL+Sex,data=cxl_F2_data))
cxl_F2_data$TDTP2_resid <- resid(lm(TDTP2~SL+Sex,data=cxl_F2_data))
cxl_F2_data$Dentary_resid <- resid(lm(Premax~SL+Sex,data=cxl_F2_data))
cxl_F2_data$Premax_resid <- resid(lm(Dentary~SL+Sex,data=cxl_F2_data))


tooth_cols = c('TVTP_resid','TDTP1_resid','TDTP2_resid','Dentary_resid','Premax_resid')
tooth_cols = c('TVTP','TDTP1','TDTP2','Dentary','Premax')
cor.test(cxl_F2_data$TVTP_resid,cxl_F2_data$Dentary_resid)
cxl_F2_cor <- cor(cxl_F2_data[tooth_cols])
cxl_F2_cor <- round(cxl_F2_cor,2)
dist <- as.dist(1-cxl_F2_cor)
hclust_order <- hclust(dist)
hclust_order 
cxl_F2_cor <- cxl_F2_cor[hclust_order$order,hclust_order$order]
cxl_F2_cor[lower.tri(cxl_F2_cor)] <- NA
cxl_F2_cor_melt <- melt(cxl_F2_cor, na.rm = TRUE)
cxl_F2_cor_melt
ggheatmap <- ggplot(cxl_F2_cor_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = .5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
Outfolder = '/home/james/Dropbox/Miller/figures/Multipop_intron4/Tooth_Counts/'
Outfolder = 'C:/Users/trahs/Dropbox/Miller/figures/Multipop_intron4/Tooth_Counts/'
ggsave(paste(Outfolder,'CxLF2s_tooth_cor_heatmap_slcorrect.pdf'))


