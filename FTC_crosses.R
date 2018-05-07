tooth_data = read.table('C:/Users/trahs/Dropbox/Miller/data/Tooth_Counts/FTC_crosses_082316.csv',header=TRUE,row.names=1,sep=',')

tooth_data = read.table('/home/james/Dropbox/Miller/data/Tooth_Counts/FTC_crosses_082316.csv',header=TRUE,row.names=1,sep=',')
ind = which(tooth_data$LVTP < tooth_data$RVTP)
max_VTPs = tooth_data$LVTP
max_VTPs[ind] <- tooth_data$RVTP[ind]  ## now contains the max's
max_VTPs 
tooth_data['max_VTP'] <- max_VTPs
plot(TVTP ~ SL, data=tooth_data)
which(tooth_data$CLUTCH == 418)
##RABS allele is in here, but HH are MD and LL are MM
tooth_data[which(tooth_data$CLUTCH == 418& tooth_data$X115_rescore == 'HL'),]$X115_rescore = NA
tooth_data[which(tooth_data$CLUTCH == 418& tooth_data$X115_rescore == 'HH'),]$X115_rescore = 'HL'
#tooth_data[which(tooth_data$CLUTCH == 511),]$TVTP = NA
#tooth_data[which(tooth_data$X115_geno == 'HH'),]$X115_geno = "HL"
backup_options <- options()
##FTC x RABS
plot(TVTP ~ SL, data=tooth_data[which(tooth_data$CLUTCH == 418),])
plot(TVTP ~ X115_rescore,data=tooth_data[which(tooth_data$CLUTCH == 418),])
points(TVTP ~ X115_rescore,data=tooth_data[which(tooth_data$CLUTCH == 418),])

##FTC x RABS
plot(TVTP ~ SL, data=tooth_data[which(tooth_data$CLUTCH == 1709),])
plot(TVTP ~ X115_rescore,data=tooth_data[which(tooth_data$CLUTCH == 1709),])
points(TVTP ~ X115_rescore,data=tooth_data[which(tooth_data$CLUTCH == 1709),])

##FTC x LITC
plot(TVTP ~ SL, data=tooth_data[which(tooth_data$CLUTCH == 511),])
plot(TVTP ~ X115_rescore,data=tooth_data[which(tooth_data$CLUTCH == 511),])
points(TVTP ~ X115_rescore,data=tooth_data[which(tooth_data$CLUTCH == 511),])

##FTC incrosses
plot(TVTP ~ SL, data=tooth_data[which(tooth_data$CLUTCH %in% c('795A','795B','795C')),])
plot(SL ~ X115_rescore, data=tooth_data[which(tooth_data$CLUTCH %in% c('795A','795B','795C')),])
plot(TVTP ~ X115_rescore,data=tooth_data[which(tooth_data$CLUTCH %in% c('795A','795B','795C')),])
summary(aov(TVTP ~ X115_rescore,data=tooth_data[which(tooth_data$CLUTCH == 511),]))

factor(FTC_marine_data$X115_rescore)
table(FTC_marine_data[!is.na(FTC_marine_data$X115_rescore),]$CLUTCH)
FTC_marine_data <- tooth_data[which(tooth_data$CLUTCH %in% c(418,511,1709)),]
FTC_marine_data['X115_rescore'] <- factor(FTC_marine_data$X115_rescore)
plot(TVTP ~ SL,data=FTC_marine_data,pch=16,col='black')
points(TVTP ~ SL,data=FTC_marine_data[FTC_marine_data$CLUTCH == 418,],pch=0)
points(TVTP ~ SL,data=FTC_marine_data[FTC_marine_data$CLUTCH == 511,],pch=1)
points(TVTP ~ SL,data=FTC_marine_data[FTC_marine_data$CLUTCH == 1709,],pch=2)

write.table(FTC_marine_data,file='/home/james/Dropbox/Miller/data/Tooth_Counts/FTC_marine_outcrosses_032218.csv',quote=FALSE,sep='\t')

points(TVTP ~ SL,data=FTC_marine_data[FTC_marine_data$X115_rescore == 'HL',],pch=16,col='blue')
points(TVTP ~ SL,data=FTC_marine_data[FTC_marine_data$X115_rescore == 'LL',],pch=16,col='red')
##Test Models
summary(lm(TVTP ~ CLUTCH,data=FTC_marine_data))
summary(lm(TVTP ~ SL,data=FTC_marine_data))
summary(lm(TVTP ~ SL+CLUTCH,data=FTC_marine_data))
summary(lm(TVTP ~ SL*CLUTCH,data=FTC_marine_data))

summary(lm(TVTP ~ SL*CLUTCH+X115_rescore,data=FTC_marine_data))
summary(lm(TVTP ~ CLUTCH+X115_rescore,data=FTC_marine_data))
AIC(lm(TVTP ~ CLUTCH,data=FTC_marine_data))
AIC(lm(TVTP ~ SL,data=FTC_marine_data))
AIC(lm(TVTP ~ SL+CLUTCH,data=FTC_marine_data))
AIC(lm(TVTP ~ SL*CLUTCH,data=FTC_marine_data))
anova(lm(TVTP ~ CLUTCH,data=FTC_marine_data),lm(TVTP ~ SL,data=FTC_marine_data),lm(TVTP ~ SL+CLUTCH,data=FTC_marine_data),lm(TVTP ~ SL*CLUTCH,data=FTC_marine_data))
anova(lm(TVTP ~ SL,data=FTC_marine_data),lm(TVTP ~ SL+CLUTCH,data=FTC_marine_data),lm(TVTP ~ SL*CLUTCH,data=FTC_marine_data))

##P~.056, but survives a one-way test
summary(aov(lm(TVTP ~ CLUTCH+X115_rescore,data=FTC_marine_data)))

options(contrasts = c('contr.sum','contr.poly'))
FTC_M_lm <- lm(TVTP ~ CLUTCH + X115_rescore, data=FTC_marine_data)
summary(FTC_M_lm)
summary(aov(FTC_M_lm))
drop1(FTC_M_lm, .~., test='F')
options(backup_options)
0.0585/2
AIC(lm(TVTP ~ X115_rescore+CLUTCH,data=FTC_marine_data))
AIC(lm(TVTP ~ X115_rescore+SL,data=FTC_marine_data))
AIC(lm(TVTP ~ X115_rescore+ SL+CLUTCH,data=FTC_marine_data))
AIC(lm(TVTP ~ X115_rescore+SL*CLUTCH,data=FTC_marine_data))

summary(lm(TVTP ~ X115_rescore,data=FTC_marine_data))
summary(lm(TVTP ~ X115_rescore+CLUTCH,data=FTC_marine_data))
summary(lm(TVTP ~ X115_rescore+SL,data=FTC_marine_data))
summary(lm(TVTP ~ X115_rescore+ SL+CLUTCH,data=FTC_marine_data))
summary(lm(TVTP ~ X115_rescore+SL*CLUTCH,data=FTC_marine_data))

ggplot(FTC_marine_data, aes(x=SL, y=TVTP,color=X115_rescore,shape=CLUTCH)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)

sl_clutch_lm <- lm(TVTP ~ CLUTCH,data=FTC_marine_data,na.action=na.exclude)
FTC_marine_data['tvtp_res'] <- resid(sl_clutch_lm,na.action=na.exclude)
FTC_marine_data <- FTC_marine_data[!is.na(FTC_marine_data$X115_rescore),]
plot(TVTP ~ X115_rescore,data=FTC_marine_data)
plot(tvtp_res ~ X115_rescore,data=FTC_marine_data)
points(tvtp_res ~ X115_rescore,data=FTC_marine_data)
summary(aov(tvtp_res ~ X115_rescore,data=FTC_marine_data))
tvtp_res_aov <- aov(tvtp_res ~ X115_rescore,data=FTC_marine_data)
TukeyHSD(tvtp_res_aov)
674.2076-673.6084
plot(tvtp_res ~ X115_rescore,data=FTC_marine_data)


#Plot the
ggplot(FTC_marine_data, aes(x=CLUTCH, y=tvtp_res,color=X115_rescore)) + 
  geom_boxplot(outlier.size=NA) +
  geom_jitter(position=position_dodge(0.8)) +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)


ggplot(FTC_marine_data, aes(x=CLUTCH, y=tvtp_res,color=X115_rescore)) + 
  geom_boxplot(outlier.size=NA) +
  geom_jitter(position=position_dodge(0.8)) +
  scale_colour_manual(values = c("HL" = "blue", 'LL' = 'darkgreen')) +
  labs(x='CLUTCH',y='Corrected Tooth Number',color='Bmp6 Genotype') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/FTC_incross_cluches.pdf',width=6,height=4.5,units='in')

ggplot(FTC_marine_data, aes(x=X115_rescore, y=tvtp_res)) + 
  geom_boxplot(outlier.size=NA) +
  geom_jitter(position=position_dodge(0.8)) +
  labs(x='Genotype',y='Corrected Tooth Number') +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size=14)) 
ggsave('/home/james/Dropbox/Miller/figures/Multipop_intron4/FTC_incross_comb.pdf',width=6,height=4.5,units='in')



FTC_marine_data[which(FTC_marine_data$X115_rescore == "LL"),]$tvtp_res
t.test(FTC_marine_data[which(FTC_marine_data$X115_rescore == "LL"),]$tvtp_res,
       FTC_marine_data[which(FTC_marine_data$X115_rescore == "HL"),]$tvtp_res,alternative="less")
wilcox.test(FTC_marine_data[which(FTC_marine_data$X115_rescore == "LL"),]$tvtp_res,
       FTC_marine_data[which(FTC_marine_data$X115_rescore == "HL"),]$tvtp_res)







##Incrosses have no effect
FTC_incrosses <- tooth_data[which(tooth_data$CLUTCH == '795A'|tooth_data$CLUTCH == '795B'|tooth_data$CLUTCH == '795C'),]
AIC(lm(TVTP ~ CLUTCH,data=FTC_incrosses))
AIC(lm(TVTP ~ SL,data=FTC_incrosses))
AIC(lm(TVTP ~ SL+CLUTCH,data=FTC_incrosses))
AIC(lm(TVTP ~ SL*CLUTCH,data=FTC_incrosses))
FTC_incrosses$tvtp_res <- resid(lm(TVTP ~ SL+CLUTCH,data=FTC_incrosses,na.action=na.exclude))
summary(aov(tvtp_res ~ X115_rescore,data=FTC_incrosses))

ggplot(FTC_incrosses, aes(x=CLUTCH, y=tvtp_res,color=X115_rescore)) + 
  geom_boxplot(outlier.size=NA) +
  geom_jitter(position=position_dodge(0.8)) +
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20)) +
  geom_smooth(method=lm,se=FALSE)
ggplot(FTC_incrosses, aes(x=SL, y=TVTP,color=X115_rescore,shape=CLUTCH)) + 
  geom_point(size=3) + 
  theme(legend.text=element_text(size=12)) + 
  theme(text = element_text(size=20))
 


table(FTC_incrosses$X115_rescore)
##This model has the best R**2, and 2nd lowest AIC (674.2 vs 673.6 for CLUTCH only)
##SL is overall positive, but the SL:CLUTCH interaction terms are negative








##Old stuff below
summary(lm(TVTP ~ CLUTCH*SL,data=tooth_data))
clutch_lm = lm(TVTP ~ CLUTCH,data=tooth_data,na.action=na.exclude)
tooth_data['corrected_TVTP'] = resid(clutch_lm,na.action=na.exclude)
plot(corrected_TVTP ~ X115_rescore,data=tooth_data)
t.test(tooth_data[which(tooth_data$X115_rescore == "LL"),]$corrected_TVTP,tooth_data[which(tooth_data$X115_rescore == "HL"|tooth_data$X115_rescore == "HH"),]$corrected_TVTP,alternative="less")


tooth_data[which(tooth_data$X115_rescore == "LL" & tooth_data$corrected_TVTP > 10),]

cl1709_data <- tooth_data[which(tooth_data$CLUTCH == 1709),]
cl1709_data[which(cl1709_data$X115_rescore == "LL"),]
t.test(cl1709_data[which(cl1709_data$X115_rescore == "LL"),]$corrected_TVTP,cl1709_data[which(cl1709_data$X115_geno == "HL"),]$corrected_TVTP)



FTC_incrosses <- tooth_data[which(tooth_data$CLUTCH == '795A'|tooth_data$CLUTCH == '795B'|tooth_data$CLUTCH == '795C'),]
plot(max_VTP ~ SL, data=FTC_incrosses)
plot(TVTP ~ SL, data=FTC_incrosses)
plot(max_VTP ~ X115_rescore, data=FTC_incrosses)
plot(TVTP ~ X115_rescore, data=FTC_incrosses)
summary(lm(TVTP ~ CLUTCH:SL, data=FTC_incrosses,na.action=na.exclude))
clutch_lm <- lm(TVTP ~ CLUTCH:SL, data=FTC_incrosses,na.action=na.exclude)
FTC_incrosses['corrected_TVTP'] = resid(clutch_lm,na.action=na.exclude)

plot(corrected_TVTP ~ X115_rescore, data=FTC_incrosses)
points(corrected_TVTP ~ X115_rescore, data=FTC_incrosses)
summary(aov(corrected_TVTP ~ X115_rescore, data=FTC_incrosses))
t.test(FTC_incrosses[which(FTC_incrosses$X115_rescore=='LL'),]$corrected_TVTP,
       FTC_incrosses[which(FTC_incrosses$X115_rescore=='HL'|FTC_incrosses$X115_rescore=='HH'),]$corrected_TVTP)


Marine_crosses <- tooth_data[which(tooth_data$CLUTCH == 511|tooth_data$CLUTCH == 1709|tooth_data$CLUTCH == 418),]
plot(TVTP ~ SL, data=Marine_crosses)
summary(lm(TVTP ~ CLUTCH+SL, data=Marine_crosses,na.action=na.exclude))
summary(lm(TVTP ~ CLUTCH:SL, data=Marine_crosses,na.action=na.exclude))
summary(lm(TVTP ~ CLUTCH*SL, data=Marine_crosses,na.action=na.exclude))
clutch_lm <- lm(TVTP ~ CLUTCH*SL, data=Marine_crosses,na.action=na.exclude)
Marine_crosses['corrected_TVTP'] = resid(clutch_lm,na.action=na.exclude)
summary(aov(corrected_TVTP ~ X115_rescore, data=Marine_crosses))
plot(TVTP ~ X115_rescore, data=Marine_crosses)
plot(corrected_TVTP ~ X115_rescore, data=Marine_crosses)
points(corrected_TVTP ~ X115_rescore, data=Marine_crosses)

Marine_crosses[which(Marine_crosses$X115_rescore=='LL'),]$corrected_TVTP
length(Marine_crosses[which(Marine_crosses$X115_rescore=='LL'),]$corrected_TVTP)
length(Marine_crosses[which(Marine_crosses$X115_rescore=='HL'|Marine_crosses$X115_rescore=='HH'),]$corrected_TVTP)
wilcox.test(Marine_crosses[which(Marine_crosses$X115_rescore=='LL'),]$corrected_TVTP,
       Marine_crosses[which(Marine_crosses$X115_rescore=='HL'|Marine_crosses$X115_rescore=='HH'),]$corrected_TVTP)
t.test(Marine_crosses[which(Marine_crosses$X115_rescore=='LL'),]$max_VTP,
       Marine_crosses[which(Marine_crosses$X115_rescore=='HL'|Marine_crosses$X115_rescore=='HH'),]$max_VTP)
t.test(Marine_crosses[which(Marine_crosses$X115_rescore=='LL'),]$TVTP,
       Marine_crosses[which(Marine_crosses$X115_rescore=='HL'|Marine_crosses$X115_rescore=='HH'),]$TVTP)


