##load plotting library
library(ggplot2)
#save options (for later, critical to run this if doing ANOVA)
backup_options <- options()
options(backup_options)

##load data
hutu_data <- read.csv('HUTU_LITC_tooth_counts.tsv',header=TRUE,sep='\t')

#two ways to plot data, default and ggplot
plot(TVTP ~ SL,data=hutu_data)
ggplot(hutu_data, aes(x=SL, y=TVTP,color=Tank)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

##look at linear models without genotype
summary(lm(TVTP ~ SL,data=hutu_data))
summary(lm(TVTP ~ Tank,data=hutu_data))
summary(lm(TVTP ~ SL + Tank,data=hutu_data))
summary(lm(TVTP ~ SL:Tank,data=hutu_data))
summary(lm(TVTP ~ SL*Tank,data=hutu_data))

##pick the model with lowest AIC
AIC(lm(TVTP ~ SL,data=hutu_data))
AIC(lm(TVTP ~ Tank,data=hutu_data))
AIC(lm(TVTP ~ SL + Tank,data=hutu_data))
AIC(lm(TVTP ~ SL:Tank,data=hutu_data))
AIC(lm(TVTP ~ SL * Tank,data=hutu_data))

##use the lowest AIC model, get residual tooth counts
tvtp_res <- resid(lm(TVTP ~ SL +Tank ,data=hutu_data,na.action=na.exclude))
hutu_data['tvtp_res'] <- tvtp_res

##plot for all markers tested
plot(tvtp_res ~ Stn.490MF,data=hutu_data)
points(tvtp_res ~ Stn.490MF,data=hutu_data)
##look at the results of a linear model, with genotype added
##Stn.490MFMF   -0.6189     1.5887  -0.390  0.69811   
##Stn.490MFMM   -4.5094     2.1453  -2.102  0.03932 * 
##the above columns test whether MF or MM is different from FF 
summary(lm(TVTP ~ SL + Tank + Stn.490MF,data=hutu_data))

##the below 4 columns do an ANOVA (don't use the aov command)
##essentially we drop each term from the model and test if the model got significantly worse
##this tests genotype as a whole, rather than comparing FF to MM
options(contrasts = c('contr.sum','contr.poly'))
tvtp_lm <- lm(TVTP ~ SL + Tank + Stn.490MF, data=hutu_data)
drop1(tvtp_lm, .~., test='F')
options(backup_options)

##more plotting for other markers
plot(tvtp_res ~ Stn.492MF,data=hutu_data)
summary(lm(TVTP ~ SL + Tank + Stn.492MF,data=hutu_data))
plot(tvtp_res ~ Stn.487,data=hutu_data)
summary(lm(TVTP ~ SL + Tank + Stn.487,data=hutu_data))

##ggplots with markers and tank, 4 variables one plot
ggplot(hutu_data, aes(x=SL, y=TVTP,color=Stn.490MF,shape=Tank)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)
ggplot(hutu_data, aes(x=SL, y=TVTP,color=Stn.492MF,shape=Tank)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)


##below here be dragons
plot(lm(TVTP ~ SL + Tank + Stn.490MF,data=hutu_data))
summary(aov(lm(TVTP ~ SL + Tank + Stn.490MF,data=hutu_data)))
summary(lm(TVTP ~ Stn.490MF,data=hutu_data))
summary(lm(TVTP ~ SL + Tank + Stn.492MF,data=hutu_data))
TukeyHSD(aov(tvtp_res ~ Stn.490MF,data=hutu_data))
summary(aov(tvtp_res ~ Stn.490MF,data=hutu_data))
summary(aov(tvtp_res ~ Stn.492MF,data=hutu_data))
summary(aov(tvtp_res ~ Stn.487,data=hutu_data))
