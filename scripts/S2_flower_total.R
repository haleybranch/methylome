####### S2 flower total
library(tidyverse)
library(lmerTest)
library(lme4) #to build the models
library(lmtest) # for likelihood ratio tests
library(ggplot2) #to graph the data
library(visreg) #to visualize residuals of the model
library(emmeans)#planned comparisons
library(multtest)#planned comparison adjustment
library(car)

rep1 <- read_csv("Data/S2_growth_rep1.csv")
rep1 <- subset(rep1,select=-c(9:17))
rep1$Rep <- 1
rep1$br_flowers <- rep1$'2flowers'
rep1$'2flowers' <- NULL
rep2 <- read_csv("Data/S2_growth_rep2.csv")
rep2 <- subset(rep2,select=-c(9:17))
rep2$Rep <- 2
rep2$br_flowers <- rep2$`2flowers`
rep2$`2flowers` <- NULL

dat <- rbind(rep1,rep2)
#s2 <- read.csv("Data/S2_IDs.csv")
s0 <- read.csv("Data/s2_s0_ids.csv")
Treatment_history <- read_csv("Data/Trthis_IDs.csv")

### merge dataframes together
dat2 <- merge(s0,Treatment_history, by="S0_ID")
dat2 <- merge(dat,dat2, by="ID")

## merge S0 and S1 treatments
dat2 <- dat2 %>%
  rowwise()%>%
  mutate(Treatment_history=paste(S0_trt,S1_trt,sep=""))
###### dead plants should have flowers as 0 

dat2 <- dat2 %>% filter(FD < 120)




#create new column with total flowers 
dat2$br_flowers <- as.numeric(dat2$br_flowers)
dat2$flower_total <- (dat2$flowers + dat2$br_flowers)

#plants that had their stem clipped, should be removed from dataset 
dat2 <- dat2 %>% filter(!is.na(flower_total))
dat2 <- dat2 %>% filter(!is.na(Diameter))


#### make all NA in flower total as 0 

## Question: Is transgeneration plasticity present?
# break up the data by region and treatment 
dat_north_wet <- filter(dat2, Region == "North", Treatment == "Wet")
dat_north_dry <- filter(dat2, Region == "North", Treatment == "Dry")
dat_south_wet <- filter(dat2, Region == "South", Treatment == "Wet")
dat_south_dry <- filter(dat2, Region == "South", Treatment == "Dry")

###################### Num of Flowers #################################################
############################################################################


################### NORTH ############################################
#### WET ######
mod_nw <- glmer.nb(flower_total~Rep + PrePeak*Treatment_history + (1|Site/S0_ID), data=dat_north_wet)

#check for model violations
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(mod_nw) # good


#anova
anova_nw <- Anova(mod_nw, type = "III") 
anova_nw
#prepeak*trthistory = 0.07 
#write.csv(anova_nw, "Results/flowertot_nw_anova.csv")


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
#write.csv(adjust.p.df, "Results/flowertot_nw_planned.csv")
#pre DW-WW = 0.02 --> DW more flowers than WW




#our predicted values from the model
predictions=predict(mod_nw,newdata=dat_north_wet, type= "response")

#the actual values
actual=dat_north_wet$flower_total

data.frame(actual=actual,predicted=predictions)%>%
  ggplot(aes(x=actual,y=predicted))+geom_point()+theme_bw()+
  geom_abline(intercept=0,slope=1)+ylim(0,30)

#makin residual plot by hand

dat_north_wet_copy <- dat_north_wet
dat_north_wet_copy$residual_vals <- predictions-dat_north_wet_copy$flower_total

ggplot(dat_north_wet_copy,aes(x=Treatment_history,y=residual_vals))+
  geom_violin()+stat_summary(fun = mean, geom="point")+facet_wrap(.~PrePeak)


vis_sla_nw<-visreg(mod_nw, xvar="Treatment_history", by="PrePeak", scale="response") 
Res_sla_nw<-vis_sla_nw$res # Extract residuals

Res_sla_nw$PrePeak <- factor(Res_sla_nw$PrePeak, level = c('Pre', 'Peak'))
Res_sla_nw$Treatment_history <- factor(Res_sla_nw$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDW", "DWW", "WDW", "WWW"))

plot_nw <- ggplot(Res_sla_nw, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  facet_wrap(.~ PrePeak) +
  scale_y_continuous(name="Flower Number")+
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
#save 8x6 


###### DRY ############

#Remove 2970 because the data does not make sense (nodes say 1.46)
dat_north_dry <- filter(dat_north_dry, ID != 2970)

mod_nd <- glmer.nb(flower_total~Rep + PrePeak*Treatment_history + (1|Site/S0_ID), data = dat_north_dry) 

#check for model violations
overdisp_fun(mod_nd) 


#anova
anova_nd <- Anova(mod_nd, type = "III") 
anova_nd
# nothing
#write.csv(anova_nd, "Results/flowertot_nd_anova.csv")

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
#nothing
#write.csv(adjust.p.df, "Results/flowertot_nd_planned.csv")

vis_sla_nd<-visreg(mod_nd, xvar="Treatment_history", by="PrePeak") 
Res_sla_nd<-vis_sla_nd$res # Extract residuals

Res_sla_nd$PrePeak <- factor(Res_sla_nd$PrePeak, level = c('Pre', 'Peak'))
Res_sla_nd$Treatment_history <- factor(Res_sla_nd$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDD", "DWD", "WDD", "WWD"))




plot_nd <- ggplot(Res_sla_nd, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
 #facet_wrap(.~ PrePeak) +
  scale_y_continuous(name="Flower Number")+
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


mod_sw <- glmer.nb(flower_total ~ Rep + PrePeak*Treatment_history + (1|Site/S0_ID), data=dat_south_wet)

#check for model violations
overdisp_fun(mod_sw)


anova_sw <- Anova(mod_sw, type = "III")
anova_sw
#nothing
#write.csv(anova_sw, "Results/flr_num_sw_anova.csv")


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
#write.csv(adjust.p.df, "Results/flr_num_sw_planned.csv")


vis_sla_sw<-visreg(mod_sw, xvar="Treatment_history", by="PrePeak") 
Res_sla_sw<-vis_sla_sw$res # Extract residuals

Res_sla_sw$PrePeak <- factor(Res_sla_sw$PrePeak, level = c('Pre', 'Peak'))
Res_sla_sw$Treatment_history <- factor(Res_sla_sw$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDW", "DWW", "WDW", "WWW"))


plot_sw <- ggplot(Res_sla_sw, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  #facet_wrap(.~ PrePeak) +
  scale_y_continuous(name="Flower Number")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_sw <-plot_sw + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_sw <-plot_sw +  #facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_sw
#save 6x4 


##### south dry ######


mod_sd <- glmer.nb(flower_total ~ Rep + PrePeak*Treatment_history + (1|Site/S0_ID), data=dat_south_dry)


#check for model violations
overdisp_fun(mod_sd) 


anova_sd <- Anova(mod_sd, type = "III")
anova_sd
#nothing
#write.csv(anova_sd, "Results/flower_num_sd_anova.csv")


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
#nothing
#write.csv(adjust.p.df, "Results/flower_num_sd_planned.csv")


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
  scale_y_continuous(name="Flower Number")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_sd <-plot_sd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_sd <-plot_sd + # facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_sd
#save 6x4






