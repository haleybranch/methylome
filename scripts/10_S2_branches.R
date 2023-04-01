###### Number of branches ######

####### S2 growth
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
rep1$br_flowers <- rep1$'2flowers'
rep1$'2flowers'<- NULL
rep2 <- read_csv("Data/S2_growth_rep2.csv")
rep2 <- subset(rep2,select=-c(9:17))
rep2$Rep <- 2
rep2$br_flowers <- rep2$`2flowers`
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

#remove dead plants from dataframe 
dat2 <- dat2 %>% filter(!is.na(Diameter))


dat2$Diameter<- as.numeric(dat2$Diameter)
dat2 <- filter(dat2, Diameter < 10)
## check sheet for dry 2621 
dat2$Rep <- as.character(dat2$Rep)

dat2$branches <- as.numeric(dat2$branches)

#remove 2300 Dry and 2681 Wet because input incorrectly
dat2 <- filter(dat2, !(ID == 2300 & Treatment == "Dry")) 
dat2 <- filter(dat2, !(ID == 2681 & Treatment == "Wet")) 

dat2$branches[is.na(dat2$branches)] <- 0

# break up the data by region and treatment 
dat_north_wet <- filter(dat2, Region == "North", Treatment == "Wet")
dat_north_dry <- filter(dat2, Region == "North", Treatment == "Dry")
dat_south_wet <- filter(dat2, Region == "South", Treatment == "Wet")
dat_south_dry <- filter(dat2, Region == "South", Treatment == "Dry")


mod_nw <- lmer(branches ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_north_wet)

#negative binomial used because count data with 0s that is overdispersed when doing poisson
mod_nw_glm <- glmer.nb(branches ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), 
                    data=dat_north_wet)


overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(mod_nw_glm)


anova_nw <- Anova(mod_nw_glm, type = "III")
anova_nw
#prepeak = 0.08
# trt history = 0.08
#prepeak * trt history = 0.04
#write.csv(anova_nw, "Results/branches_nw_anova.csv")



###
emm1 = emmeans(mod_nw_glm, specs = ~ PrePeak*Treatment_history)
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
#write.csv(adjust.p.df, "Results/branch_nw_planned.csv")
#pre DW-WW < 0.001 --> WW fewer than DW
#pre DD-WW = 0.001 --> WW fewer than DD
#pre-peak WW = 0.002 --> peak more than pre 
#peak DD-WD = 0.05 --> WD fewer than DD

vis_sla_nw<-visreg(mod_nw, xvar="Treatment_history", by="PrePeak") 
Res_sla_nw<-vis_sla_nw$res # Extract residuals

Res_sla_nw$PrePeak <- factor(Res_sla_nw$PrePeak, level = c('Pre', 'Peak'))
Res_sla_nw$Treatment_history <- factor(Res_sla_nw$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDW", "DWW", "WDW", "WWW"))

summary(mod_nw)
plot_nw <- ggplot(Res_sla_nw, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="skyblue3", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  facet_wrap(.~ PrePeak) +
  scale_y_continuous(name="Number of Branches")+
  #scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  #scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_nw <-plot_nw + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_nw <-plot_nw + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_nw
#save 6x4 

##### north dry ######


mod_nd <- glmer.nb(branches ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_north_dry)



overdisp_fun(mod_nd)


anova_nd <- Anova(mod_nd, type = "III")
anova_nd
#nothing
#write.csv(anova_nd, "Results/branches_nd_anova.csv")


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
#no differences
#write.csv(adjust.p.df, "Results/branches_nd_planned.csv")

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
  scale_y_continuous(name="Number of Branches")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_nd <-plot_nd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_nd <-plot_nd + #facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_nd
#save 6x4 



################# SOUTH ##################


mod_sw <- glmer.nb(branches ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_south_wet)

overdisp_fun(mod_sw)


anova_sw <- Anova(mod_sw, type = "III")
anova_sw
#nothing
#write.csv(anova_sw, "Results/branches_sw_anova.csv")


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
#nothing
#write.csv(adjust.p.df, "Results/branches_sw_planned.csv")

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
  scale_y_continuous(name="Number of Branches")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
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
#save 6x4 


##### south dry ######


mod_sd <- glmer.nb(branches ~ Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_south_dry)

overdisp_fun(mod_sd)


anova_sd <- Anova(mod_sd, type = "III")
anova_sd
#nothing
#write.csv(anova_sd, "Results/branches_sd_anova.csv")


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
#write.csv(adjust.p.df, "Results/branches_sd_planned.csv")


vis_sla_sd<-visreg(mod_sd, xvar="Treatment_history", by="PrePeak") 
Res_sla_sd<-vis_sla_sd$res # Extract residuals

Res_sla_sd$PrePeak <- factor(Res_sla_sd$PrePeak, level = c('Pre', 'Peak'))
Res_sla_sd$Treatment_history <- factor(Res_sla_sd$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDD", "DWD", "WDD", "WWD"))


plot_sd <- ggplot(Res_sla_sd, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
#  facet_wrap(.~ PrePeak) +
  scale_y_continuous(name="Number of Branches")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_sd <-plot_sd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_sd <-plot_sd + #facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_sd
#save 8x6 
