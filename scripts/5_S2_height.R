####### plant height ########
library(tidyverse)
library(lmerTest)
library(lme4) #to build the models
library(ggplot2) #to graph the data
library(visreg) #to visualize residuals of the model
library(emmeans)#planned comparisons
library(multtest)#planned comparison adjustment
library(effectsize)


rep1 <- read_csv("Data/S2_growth_rep1.csv")
rep1 <- subset(rep1,select=-c(9:17))
rep1$Rep <- 1
rep1$'2flowers' <- NULL

rep2 <- read_csv("Data/S2_growth_rep2.csv")
rep2 <- subset(rep2,select=-c(9:17))
rep2$Rep <- 2
rep2$`2flowers` <- NULL

dat <- rbind(rep1,rep2)
s0 <- read.csv("Data/s2_s0_ids.csv")
Treatment_history <- read_csv("Data/Trthis_IDs.csv")

### merge dataframes together
dat2 <- merge(s0,Treatment_history, by="S0_ID")
dat2 <- merge(dat,dat2, by="ID")

## merge S0 and S1 treatments
dat2 <- dat2 %>%
  rowwise()%>%
  mutate(Treatment_history=paste(S0_trt,S1_trt,sep=""))

dat2 <- dat2 %>% filter(!is.na(Height))
dat2$Height<- as.numeric(dat2$Height)
dat2$Rep <- as.character(dat2$Rep)

# break up the data by region and treatment 
dat_north_wet <- filter(dat2, Region == "North", Treatment == "Wet")
dat_north_dry <- filter(dat2, Region == "North", Treatment == "Dry")
dat_south_wet <- filter(dat2, Region == "South", Treatment == "Wet")
dat_south_dry <- filter(dat2, Region == "South", Treatment == "Dry")



############## NORTH WET ##########3

mod_nw <- lmer(Height ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_north_wet)

#check for model violations
plot(mod_nw, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_nw,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_nw) # qqplot

hist(resid(mod_nw)) #histogram

plot(mod_nw, rstudent(.) ~ hatvalues(.)) #residuals vs leverage


anova_nw <- anova(mod_nw)
anova_nw
#treatment history = 0.09
#write.csv(anova_nw, "Results/height_nw_anova.csv")


#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_nw <- omega_squared(mod_nw)
effsize_nw
#TRTHIS < 0
#write.csv(effsize_nw, "Results/height_nw_effectsize.csv")


###
emm1 = emmeans(mod_nw, specs = ~ PrePeak*Treatment_history)
emm1


peak_DD = c(1, 0, 0, 0, 0, 0, 0, 0)
pre_DD = c(0, 1, 0, 0, 0, 0, 0, 0)
peak_DW = c(0, 0, 1, 0, 0, 0, 0, 0)
pre_DW = c(0, 0, 0, 1, 0, 0, 0, 0)
peak_WD = c(0, 0, 0, 0, 1, 0, 0, 0)
pre_WD = c(0, 0, 0, 0, 0, 1, 0, 0)
peak_WW = c(0, 0, 0, 0, 0, 0, 1, 0)
pre_WW = c(0, 0, 0, 0, 0, 0, 0, 1)



contrasts <- list("pre_DD-DW" = pre_DD - pre_DW, "pre_DD-WD" = pre_DD - pre_WD,
                  "pre_DD-WW" = pre_DD - pre_WW, 
                  "pre_DW-WW" = pre_DW - pre_WW, "pre_WD-WW" = pre_WD - pre_WW,
                  "peak_DD-DW" = peak_DD - peak_DW, "peak_DD-WD" = peak_DD - peak_WD,
                  "peak_DD-WW" = peak_DD - peak_WW, 
                  "peak_DW-WW" = peak_DW - peak_WW, "peak_WD-WW" = peak_WD - peak_WW,
                  
                  "pre-peak_DD" = pre_DD - peak_DD, "pre-peak_DW" = pre_DW - peak_DW,
                  "pre-peak_WD" = pre_WD - peak_WD, "pre-peak_WW" = pre_WW - peak_WW)


#grab the raw pvalues without a posthoc adjustment
contrast.amax <- contrast(emm1, method = contrasts, adjust="none") 
summ.contrast <- summary(contrast.amax)

#use a pvalue correction for the multiple tests 
adjust.p <- mt.rawp2adjp(summ.contrast$p.value, proc=c("BH"), alpha = 0.05, na.rm = FALSE)
#adjust.p <- as.data.frame(adjust.p$adjp)
#adjust.p <- filter(adjust.p, rawp > 0.0001)
#realign the adjusted pvalues with the contrast names
adjust.p<-as.data.frame(adjust.p$adjp)
adjust.p<-distinct(adjust.p)
adjust.p.df<- adjust.p %>%
  rename(p.value=rawp)%>%
  inner_join(summ.contrast)
adjust.p.df
#write.csv(adjust.p.df, "Results/height_nw_planned.csv")
#nothing

vis_sla_nw<-visreg(mod_nw, xvar="Treatment_history", by="PrePeak") 
Res_sla_nw<-vis_sla_nw$res # Extract residuals

Res_sla_nw$PrePeak <- factor(Res_sla_nw$PrePeak, level = c('Pre', 'Peak'))
Res_sla_nw$Treatment_history <- factor(Res_sla_nw$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDW", "DWW", "WDW", "WWW"))

plot_nw <- ggplot(Res_sla_nw, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="skyblue3", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  facet_wrap(.~ PrePeak) +
  ylab("Plant Height (cm)")+
  ylim(0,100)+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_nw <-plot_nw + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_nw <-plot_nw +  facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_nw
#save 6x4


##### north dry ######


mod_nd <- lmer(Height ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_north_dry)

#check for model violations
plot(mod_nd, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_nd,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_nd) # qqplot


hist(resid(mod_nd)) #histogram

plot(mod_nd, rstudent(.) ~ hatvalues(.)) #residuals vs leverage


anova_nd <- anova(mod_nd)
anova_nd
#nothing
#write.csv(anova_nd, "Results/height_nd_anova.csv")

#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_nd <- omega_squared(mod_nd)
effsize_nd
#TRTHIS < 0
#write.csv(effsize_nd, "Results/height_nd_effectsize.csv")

###
emm1 = emmeans(mod_nd, specs = ~ PrePeak*Treatment_history)
emm1


peak_DD = c(1, 0, 0, 0, 0, 0, 0, 0)
pre_DD = c(0, 1, 0, 0, 0, 0, 0, 0)
peak_DW = c(0, 0, 1, 0, 0, 0, 0, 0)
pre_DW = c(0, 0, 0, 1, 0, 0, 0, 0)
peak_WD = c(0, 0, 0, 0, 1, 0, 0, 0)
pre_WD = c(0, 0, 0, 0, 0, 1, 0, 0)
peak_WW = c(0, 0, 0, 0, 0, 0, 1, 0)
pre_WW = c(0, 0, 0, 0, 0, 0, 0, 1)



contrasts <- list("pre_DD-DW" = pre_DD - pre_DW, "pre_DD-WD" = pre_DD - pre_WD,
                  "pre_DD-WW" = pre_DD - pre_WW, 
                  "pre_DW-WW" = pre_DW - pre_WW, "pre_WD-WW" = pre_WD - pre_WW,
                  "peak_DD-DW" = peak_DD - peak_DW, "peak_DD-WD" = peak_DD - peak_WD,
                  "peak_DD-WW" = peak_DD - peak_WW, 
                  "peak_DW-WW" = peak_DW - peak_WW, "peak_WD-WW" = peak_WD - peak_WW,
                  
                  "pre-peak_DD" = pre_DD - peak_DD, "pre-peak_DW" = pre_DW - peak_DW,
                  "pre-peak_WD" = pre_WD - peak_WD, "pre-peak_WW" = pre_WW - peak_WW)


#grab the raw pvalues without a posthoc adjustment
contrast.amax <- contrast(emm1, method = contrasts, adjust="none") 
summ.contrast <- summary(contrast.amax)

#use a pvalue correction for the multiple tests 
adjust.p <- mt.rawp2adjp(summ.contrast$p.value, proc=c("BH"), alpha = 0.05, na.rm = FALSE)
#adjust.p <- as.data.frame(adjust.p$adjp)
#adjust.p <- filter(adjust.p, rawp > 0.0001)
#realign the adjusted pvalues with the contrast names
adjust.p<-as.data.frame(adjust.p$adjp)
adjust.p<-distinct(adjust.p)
adjust.p.df<- adjust.p %>%
  rename(p.value=rawp)%>%
  inner_join(summ.contrast)
adjust.p.df
#write.csv(adjust.p.df, "Results/height_nd_planned.csv")
#nothing

vis_sla_nd<-visreg(mod_nd, xvar="Treatment_history", by="PrePeak") 
Res_sla_nd<-vis_sla_nd$res # Extract residuals

Res_sla_nd$PrePeak <- factor(Res_sla_nd$PrePeak, level = c('Pre', 'Peak'))
Res_sla_nd$Treatment_history <- factor(Res_sla_nd$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDD", "DWD", "WDD", "WWD"))


plot_nd <- ggplot(Res_sla_nd, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  facet_wrap(.~ PrePeak) +
  ylab("Plant Height (cm)")+
  ylim(0,100)+  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_nd <-plot_nd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_nd <-plot_nd + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_nd
#save 6x4 




################# SOUTH ##################


mod_sw <- lmer(Height ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_south_wet)

#check for model violations
plot(mod_sw, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_sw,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_sw) # qqplot

hist(resid(mod_sw)) #histogram

plot(mod_sw, rstudent(.) ~ hatvalues(.)) #residuals vs leverage


anova_sw <- anova(mod_sw)
anova_sw
#PrePeak*TrtHistory = 0.008
#treatment history = 0.005
#prepeak = 0.05
#write.csv(anova_sw, "Results/height_sw_anova.csv")


#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_sw <- omega_squared(mod_sw)
effsize_sw
#TRTHIS = 0.03
#prepeak*trthis = 0.02
#prepeak = 0.15
#write.csv(effsize_sw, "Results/height_sw_effectsize.csv")


###
emm1 = emmeans(mod_sw, specs = ~ PrePeak*Treatment_history)
emm1


peak_DD = c(1, 0, 0, 0, 0, 0, 0, 0)
pre_DD = c(0, 1, 0, 0, 0, 0, 0, 0)
peak_DW = c(0, 0, 1, 0, 0, 0, 0, 0)
pre_DW = c(0, 0, 0, 1, 0, 0, 0, 0)
peak_WD = c(0, 0, 0, 0, 1, 0, 0, 0)
pre_WD = c(0, 0, 0, 0, 0, 1, 0, 0)
peak_WW = c(0, 0, 0, 0, 0, 0, 1, 0)
pre_WW = c(0, 0, 0, 0, 0, 0, 0, 1)



contrasts <- list("pre_DD-DW" = pre_DD - pre_DW, "pre_DD-WD" = pre_DD - pre_WD,
                  "pre_DD-WW" = pre_DD - pre_WW, 
                  "pre_DW-WW" = pre_DW - pre_WW, "pre_WD-WW" = pre_WD - pre_WW,
                  "peak_DD-DW" = peak_DD - peak_DW, "peak_DD-WD" = peak_DD - peak_WD,
                  "peak_DD-WW" = peak_DD - peak_WW, 
                  "peak_DW-WW" = peak_DW - peak_WW, "peak_WD-WW" = peak_WD - peak_WW,
                  
                  "pre-peak_DD" = pre_DD - peak_DD, "pre-peak_DW" = pre_DW - peak_DW,
                  "pre-peak_WD" = pre_WD - peak_WD, "pre-peak_WW" = pre_WW - peak_WW)


#grab the raw pvalues without a posthoc adjustment
contrast.amax <- contrast(emm1, method = contrasts, adjust="none") 
summ.contrast <- summary(contrast.amax)

#use a pvalue correction for the multiple tests 
adjust.p <- mt.rawp2adjp(summ.contrast$p.value, proc=c("BH"), alpha = 0.05, na.rm = FALSE)
#adjust.p <- as.data.frame(adjust.p$adjp)
#adjust.p <- filter(adjust.p, rawp > 0.0001)
#realign the adjusted pvalues with the contrast names
adjust.p<-as.data.frame(adjust.p$adjp)
adjust.p<-distinct(adjust.p)
adjust.p.df<- adjust.p %>%
  rename(p.value=rawp)%>%
  inner_join(summ.contrast)
adjust.p.df
#write.csv(adjust.p.df, "Results/height_sw_planned.csv")
#peak DD-WW < 0.001 -->
#peak DW-WW = 0.008 -->
#pre-peak DD = 0.009 -->
#peak DD-WD = 0.009 -->
#pre-peak DW = 0.07 --> 

vis_sla_sw<-visreg(mod_sw, xvar="Treatment_history", by="PrePeak") 
Res_sla_sw<-vis_sla_sw$res # Extract residuals

Res_sla_sw$PrePeak <- factor(Res_sla_sw$PrePeak, level = c('Pre', 'Peak'))
Res_sla_sw$Treatment_history <- factor(Res_sla_sw$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDW", "DWW", "WDW", "WWW"))


plot_sw <- ggplot(Res_sla_sw, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="skyblue3", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  facet_wrap(.~ PrePeak) +
  ylab("Plant Height (cm)")+
  ylim(0,100)+  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_sw <-plot_sw + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_sw <-plot_sw + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_sw
#save 8x6 


##### south dry ######


mod_sd <- lmer(Height ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_south_dry)

#check for model violations
plot(mod_sd, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_sd,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_sd) # qqplot

hist(resid(mod_sd)) #histogram

plot(mod_sd, rstudent(.) ~ hatvalues(.)) #residuals vs leverage


anova_sd <- anova(mod_sd)
anova_sd
#prepeak*treatment history < 0.0001
#treatment history = 0.08
#write.csv(anova_sd, "Results/height_sd_anova.csv")



#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_sd <- omega_squared(mod_sd)
effsize_sd
#TRTHIS < 0
# Prepeak*trt his = 0.05
#write.csv(effsize_sd, "Results/height_sd_effectsize.csv")


###
emm1 = emmeans(mod_sd, specs = ~ PrePeak*Treatment_history)
emm1


peak_DD = c(1, 0, 0, 0, 0, 0, 0, 0)
pre_DD = c(0, 1, 0, 0, 0, 0, 0, 0)
peak_DW = c(0, 0, 1, 0, 0, 0, 0, 0)
pre_DW = c(0, 0, 0, 1, 0, 0, 0, 0)
peak_WD = c(0, 0, 0, 0, 1, 0, 0, 0)
pre_WD = c(0, 0, 0, 0, 0, 1, 0, 0)
peak_WW = c(0, 0, 0, 0, 0, 0, 1, 0)
pre_WW = c(0, 0, 0, 0, 0, 0, 0, 1)



contrasts <- list("pre_DD-DW" = pre_DD - pre_DW, "pre_DD-WD" = pre_DD - pre_WD,
                  "pre_DD-WW" = pre_DD - pre_WW, 
                  "pre_DW-WW" = pre_DW - pre_WW, "pre_WD-WW" = pre_WD - pre_WW,
                  "peak_DD-DW" = peak_DD - peak_DW, "peak_DD-WD" = peak_DD - peak_WD,
                  "peak_DD-WW" = peak_DD - peak_WW, 
                  "peak_DW-WW" = peak_DW - peak_WW, "peak_WD-WW" = peak_WD - peak_WW,
                  
                  "pre-peak_DD" = pre_DD - peak_DD, "pre-peak_DW" = pre_DW - peak_DW,
                  "pre-peak_WD" = pre_WD - peak_WD, "pre-peak_WW" = pre_WW - peak_WW)


#grab the raw pvalues without a posthoc adjustment
contrast.amax <- contrast(emm1, method = contrasts, adjust="none") 
summ.contrast <- summary(contrast.amax)

#use a pvalue correction for the multiple tests 
adjust.p <- mt.rawp2adjp(summ.contrast$p.value, proc=c("BH"), alpha = 0.05, na.rm = FALSE)
#adjust.p <- as.data.frame(adjust.p$adjp)
#adjust.p <- filter(adjust.p, rawp > 0.0001)
#realign the adjusted pvalues with the contrast names
adjust.p<-as.data.frame(adjust.p$adjp)
adjust.p<-distinct(adjust.p)
adjust.p.df<- adjust.p %>%
  rename(p.value=rawp)%>%
  inner_join(summ.contrast)
adjust.p.df
#write.csv(adjust.p.df, "Results/height_sd_planned.csv")
#peak DW-WW -->
#pre-peak DW --> 

vis_sla_sd<-visreg(mod_sd, xvar="Treatment_history", by="PrePeak") 
Res_sla_sd<-vis_sla_sd$res # Extract residuals

Res_sla_sd$PrePeak <- factor(Res_sla_sd$PrePeak, level = c('Pre', 'Peak'))
Res_sla_sd$Treatment_history <- factor(Res_sla_sd$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDD", "DWD", "WDD", "WWD"))


plot_sd <- ggplot(Res_sla_sd, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  facet_wrap(.~ PrePeak) +
  ylab("Plant Height (cm)")+
  ylim(0,100)+  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_sd <-plot_sd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_sd <-plot_sd + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_sd
#save 8x6 
