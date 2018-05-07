library(ggplot2)
library(reshape2)
backup_options <- options()

Folder = '/home/james/Dropbox/Miller/data/Mutants/Bmp6/'
Folder = 'C:/Users/trahs/Dropbox/Miller/data/Mutants/Bmp6/'
Dropbox = 'C:/Users/trahs/Dropbox/'

bmp6_oral <- read.csv(paste(Folder,'Bmp6_exonTALENs_all_072315.csv',sep=''))
bmp6_oral <- bmp6_oral[!is.na(bmp6_oral$den),]               
bmp6_oral <- bmp6_oral[!is.na(bmp6_oral$pm),]   
bmp6_oral <- bmp6_oral[bmp6_oral$pm < 60,]   
bmp6_oral <- bmp6_oral[!is.na(bmp6_oral$genotype),]   
bmp6_oral$genotype <-  factor(bmp6_oral$genotype,levels=c('WW','WD'))


mean(bmp6_oral$pm)
sd(bmp6_oral[bmp6_oral$genotype == 'WW',]$pm)

ggplot(bmp6_oral, aes(x=SL, y=TVTP,color=genotype)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(bmp6_oral, aes(x=SL, y=cb1.mm,color=genotype)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(bmp6_oral, aes(x=SL, y=cb2.mm,color=genotype)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(bmp6_oral, aes(x=SL, y=cb3.mm,color=genotype)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(bmp6_oral, aes(x=SL, y=cb4.mm,color=genotype)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

ggplot(bmp6_oral, aes(x=SL, y=cb5.mm,color=genotype)) + 
  scale_colour_manual(values = c("WW" = "blue", 'WD' = 'purple', "DD" = "red")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method='lm',se=FALSE)

summary(lm(TVTP ~ genotype + SL + cb5.mm,data=bmp6_oral))

summary(lm(cb1.mm ~ genotype + SL,data=bmp6_oral))
summary(lm(cb2.mm ~ genotype + SL,data=bmp6_oral))
summary(lm(cb3.mm ~ genotype + SL,data=bmp6_oral))
summary(lm(cb4.mm ~ genotype + SL,data=bmp6_oral))
summary(lm(cb5.mm ~ genotype + SL,data=bmp6_oral))
0.0615/2
summary(lm(pm ~ genotype+SL,data=bmp6_oral))
tooth_lm <- lm(pm ~ genotype+SL,data=bmp6_oral)
tooth_predict <- predict(tooth_lm,interval="confidence")
colnames(tooth_predict) <- c('pm_fit','pm_low','pm_up')
bmp6_oral<-cbind(bmp6_oral,tooth_predict)
ggplot(bmp6_oral, aes(x=SL, y=pm,color=genotype)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=pm_fit),size=2)
ggsave(paste(Dropbox,'Miller/figures/Mutants/Bmp6/Bmp6_oral_SL_pm.pdf',sep=""))

bmp6_oral$pm_resid <- resid(lm(pm ~ SL,data=bmp6_oral))
ggplot(bmp6_oral,aes(x=genotype,y=pm_resid,color=genotype)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Premaxilla') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Bmp6 \nGenotype",override.aes = list(size=3,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Bmp6/Bmp6_oral_box_pm_resid.pdf',sep=""))
##excluding outlier still is signif by a one-tailed test
t.test(bmp6_oral[bmp6_oral$genotype == 'WW',]$pm_resid,
       bmp6_oral[bmp6_oral$genotype == 'WD',]$pm_resid,alternative='greater')

summary(lm(den ~ genotype+SL,data=bmp6_oral))
tooth_lm <- lm(den ~ genotype+SL,data=bmp6_oral)
tooth_predict <- predict(tooth_lm,interval="confidence")
colnames(tooth_predict) <- c('den_fit','den_low','den_up')
bmp6_oral<-cbind(bmp6_oral,tooth_predict)
ggplot(bmp6_oral, aes(x=SL, y=den,color=genotype)) + 
  scale_colour_manual(values = c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  geom_point(size=3) + 
  theme_bw() +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_line(aes(y=den_fit),size=2)
ggsave(paste(Dropbox,'Miller/figures/Mutants/Bmp6/Bmp6_oral_SL_den.pdf',sep=""))


bmp6_oral$den_resid <- resid(lm(den ~ SL,data=bmp6_oral))
ggplot(bmp6_oral,aes(x=genotype,y=den_resid,color=genotype)) +
  #geom_violin(outlier.shape=NA) +
  geom_boxplot(outlier.shape=NA,size=1.5) +
  geom_jitter(size=3) +
  scale_colour_manual(values =  c("DD" = "coral", 'WD' = 'magenta', "WW" = "royalblue")) +
  theme_bw() +
  labs(y='Corrected Dentary') +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=30)) +
  guides(color = guide_legend(override.aes = list(size=3,linetype=0))) + 
  guides(color=guide_legend(title="Bmp6 \nGenotype",override.aes = list(size=3,linetype=0)))
ggsave(paste(Dropbox,'Miller/figures/Mutants/Bmp6/Bmp6_oral_box_den_resid.pdf',sep=""))

t.test(bmp6_oral[bmp6_oral$genotype == 'WW',]$den_resid,
       bmp6_oral[bmp6_oral$genotype == 'WD',]$den_resid,alternative='greater')

summary(lm(den ~ genotype + SL,data=bmp6_oral))
summary(lm(pm ~ genotype + SL,data=bmp6_oral))
