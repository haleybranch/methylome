library(tidyverse)
library(lme4)
library(lmerTest)


seed <- read.csv("S0_S1_S2_seed_aug16.csv")
seed <- na.omit(seed)
seed$Weight_1seed <- as.numeric(seed$Weight_1seed)
seed$Num_seeds <- as.numeric(seed$Num_seeds)

seed$ID <- as.character(seed$ID)
seed <- seed %>%
  rowwise()%>%
  mutate(TRTHIS=paste(Parental.ignore.,Treatment,sep=""))

seed <- seed %>% filter(Weight_1seed < 20)

boxplot(seed$Weight_1seed ~ seed$TRTHIS)
seed$S0_ID <- NULL
names(seed)[names(seed) == 'Timeline..ignore.'] <- 'S0_ID'

ids <- read.csv("Data/S2_s0_ids.csv")
new <- inner_join(seed,ids, by="S0_ID")


mod1 <- lmer(Weight_1seed ~ Treatment*PrePeak*Region*Parental.ignore. + (1|Site/S0_ID), data=new)
anova(mod1)

pre_FT_d <-visreg(mod1, xvar="Parental.ignore.", by="Region", cond=list(Treatment="dry", PrePeak = "Pre")) #set up visreg for Drought
peak_FT_d <-visreg(mod1, xvar="Parental.ignore.", by="Region", cond=list(Treatment="dry", PrePeak = "Peak")) #set up visreg for Drought

pre_FT_w<-visreg(mod1, xvar="Parental.ignore.", by="Region", cond=list(Treatment="wet", PrePeak = "Pre")) #set up visreg for Drought
peak_FT_w<-visreg(mod1, xvar="Parental.ignore.", by="Region", cond=list(Treatment="wet", PrePeak = "Peak")) #set up visreg for Drought

Res_FT_pre_d<-pre_FT_d$res; Res_FT_pre_w<-pre_FT_w$res # Extract residuals
Res_FT_peak_d<-peak_FT_d$res; Res_FT_peak_w<-peak_FT_w$res # Extract residuals
Res_FT_peak <- rbind(Res_FT_pre_d,Res_FT_pre_w,Res_FT_peak_d,Res_FT_peak_w)


# plot
plot <-ggplot(Res_FT_peak, aes(Parental.ignore., y=visregRes, colour=Region))+
  #geom_boxplot(aes(colour=Treatment)) +
  geom_violin(aes(fill=Treatment), position = position_dodge(width = 1))+
  stat_summary(aes(group=Treatment, colour=Treatment), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(group=Treatment, colour=Treatment),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  facet_grid(PrePeak ~ Region) +
  scale_y_continuous(name="Seed Weight")+
  scale_color_manual(values= c("dry"="#FFA100", "wet"="skyblue3"))+
  scale_fill_manual(values= c("dry"="#FFA100", "wet"="skyblue3"))+
  theme_classic()
plot <-plot + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot <-plot + #facet_wrap(.~Region) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot




emm1 = emmeans(mod1, specs = ~ Treatment*Region*PrePeak*TRTHIS, adjust="Tukey")
emm1


DNPe = c(1, 0, 0, 0, 0, 0, 0, 0)
WNPe = c(0, 1, 0, 0, 0, 0, 0, 0)
DSPe = c(0, 0, 1, 0, 0, 0, 0, 0)
WSPe = c(0, 0, 0, 1, 0, 0, 0, 0)
DNPr = c(0, 0, 0, 0, 1, 0, 0, 0)
WNPr = c(0, 0, 0, 0, 0, 1, 0, 0)
DSPr = c(0, 0, 0, 0, 0, 0, 1, 0)
WSPr = c(0, 0, 0, 0, 0, 0, 0, 1)

contrasts <- list("DNPe-WNPe" = DNPe - WNPe, "DSPe-WSPe" = DSPe - WSPe,
                  "DNPr-WNPr" = DNPr - WNPr, "DSPr-WSPr" = DSPr - WSPr,
                  "DNPe-DNPr" = DNPe - DNPr, "WNPe-WNPr" = WNPe - WNPr,
                  "DSPe-DSPr" = DSPe - DSPr, "WSPe-WSPr" = WSPe - WSPr,
                  "DNPe-DSPe" = DNPe - DSPe, "WNPe-WSPe" = WNPe - WSPe,
                  "DNPr-DSPr" = DNPr - DSPr, "WNPr-WSPr" = WNPr - WSPr)

contrast(emm1, method = contrasts, adjust="none")


#filter region 
south_seed_peak_dd <- filter(new, Region=="South", PrePeak == "Peak", TRTHIS=="Ddry")
south_seed_peak_dw <- filter(new, Region=="South", PrePeak == "Peak", TRTHIS=="Dwet")
south_seed_peak_wd <- filter(new, Region=="South", PrePeak == "Peak", TRTHIS=="Wdry")
south_seed_peak_ww <- filter(new, Region=="South", PrePeak == "Peak", TRTHIS=="Wwet")

south_seed_pre_dd <- filter(new, Region=="South", PrePeak == "Pre", TRTHIS=="Ddry")
south_seed_pre_dw <- filter(new, Region=="South", PrePeak == "Pre", TRTHIS=="Dwet")
south_seed_pre_wd <- filter(new, Region=="South", PrePeak == "Pre", TRTHIS=="Wdry")
south_seed_pre_ww <- filter(new, Region=="South", PrePeak == "Pre", TRTHIS=="Wwet")

north_seed_peak_dd <- filter(new, Region=="North", PrePeak == "Peak", TRTHIS=="Ddry")
north_seed_peak_dw <- filter(new, Region=="North", PrePeak == "Peak", TRTHIS=="Dwet")
north_seed_peak_wd <- filter(new, Region=="North", PrePeak == "Peak", TRTHIS=="Wdry")
north_seed_peak_ww <- filter(new, Region=="North", PrePeak == "Peak", TRTHIS=="Wwet")

north_seed_pre_dd <- filter(new, Region=="North", PrePeak == "Pre", TRTHIS=="Ddry")
north_seed_pre_dw <- filter(new, Region=="North", PrePeak == "Pre", TRTHIS=="Dwet")
north_seed_pre_wd <- filter(new, Region=="North", PrePeak == "Pre", TRTHIS=="Wdry")
north_seed_pre_ww <- filter(new, Region=="North", PrePeak == "Pre", TRTHIS=="Wwet")


library(ggplot2)
Site_Labs<-c("2010" = "Pre", )
ggplot(seed, aes(trthis, y=))



plot<- ggplot(new, aes(TRTHIS, y=ratio, colour=PrePeak))+
  geom_boxplot(aes(colour=Region))+
  xlab("Treatment History") +
  facet_wrap(.~ Region) +
  theme_classic()
plot <-plot + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot <-plot + facet_wrap(.~Region) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot

mean(seed$Weight_1seed)

new$avg_seed <- new$Weight_total..mg./0.02366806

new$Num_seeds <- as.numeric(new$Num_seeds)

plot(new$Num_seeds ~ new$avg_seed)

new$ratio <- (new$Num_seeds)/(new$avg_seed)

plot(new$Num_seeds~new$avg_seed)
