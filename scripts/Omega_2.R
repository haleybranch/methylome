### Comparing effect sizes for year and treatment history #####

###
library(dplyr)

north_wet <- read.csv("Data/ES_north_wet.csv")
north_dry <- read.csv("Data/ES_north_dry.csv")
south_wet <- read.csv("Data/ES_south_wet.csv")
south_dry <- read.csv("Data/ES_south_dry.csv")

### North Wet ####

north_wet %>% 
  group_by(Factor) %>%
  summarise(
    count = n(),
    mean = mean(ES),
    sd = sd(ES),
  )

differences <- with(north_wet, ES[Factor == "Year"] - ES[Factor == "TRTHIS"])
shapiro.test(differences)

north_wet$ES <- north_wet$ES + 0.0001
north_wet$Trait <- factor(north_wet$Trait, level = c("SLA", "Leaf #", "Diameter", "Height", "Biomass","Germ","Flowering","Flower #", "Root"), labels = c("SLA", "Leaf #", "Diameter", "Height", "Biomass","Germ","Flowering","Flower #", "Root"))

n_wet <- ggplot(north_wet, aes(Trait, y=ES, fill = Factor))+
  geom_bar(stat="identity", position = "dodge")+
  ggtitle("North Wet")+
  ylim(0, 0.19) +
  scale_color_manual(values= c("Year"="#FFA100", "TRTHIS"="skyblue3"))+
  scale_fill_manual(values= c("Year"="#FFA100", "TRTHIS"="skyblue3"))+
  theme_classic() +
  theme(legend.position = "none")
n_wet <-n_wet + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
  title = element_text(color="black", size=30))
n_wet

# no difference

#### North Dry #####

north_dry %>% 
  group_by(Factor) %>%
  summarise(
    count = n(),
    mean = mean(ES),
    sd = sd(ES),
  )

north_dry$ES <- north_dry$ES + 0.0001
north_dry$Trait <- factor(north_dry$Trait, level = c("SLA", "Leaf #", "Diameter", "Height", "Biomass","Germ","Flowering","Flower #", "Root"), labels = c("SLA", "Leaf #", "Diameter", "Height", "Biomass","Germ","Flowering","Flower #", "Root"))

n_dry <- ggplot(north_dry, aes(Trait, y=ES, fill = Factor))+
  geom_bar(stat="identity", position = "dodge")+
  ggtitle("North Dry")+
  ylim(0, 0.19) +
  scale_color_manual(values= c("Year"="#FFA100", "TRTHIS"="skyblue3"))+
  scale_fill_manual(values= c("Year"="#FFA100", "TRTHIS"="skyblue3"))+
  theme_classic() +
  theme(legend.position = "none")
n_dry <-n_dry + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
  title = element_text(color="black", size=30))
n_dry


north_dry$ES <- north_dry$ES + 0.0001

south_wet$ES <- south_wet$ES + 0.0001
library(ggplot2)

south_wet$Trait <- factor(south_wet$Trait, level = c("SLA", "Leaf #", "Diameter", "Height", "Biomass","Germ","Flowering","Flower #", "Root"), labels = c("SLA", "Leaf #", "Diameter", "Height", "Biomass","Germ","Flowering","Flower #", "Root"))


s_w <- ggplot(south_wet, aes(Trait, y=ES, fill = Year))+
  geom_bar(stat="identity", position = "dodge")+
  ggtitle("South Wet")+
  ylim(0, 0.19) +
  scale_color_manual(values= c("Year"="#FFA100", "TRTHIS"="skyblue3"))+
  scale_fill_manual(values= c("Year"="#FFA100", "TRTHIS"="skyblue3"))+
  theme_classic() +
  theme(legend.position = "none")
s_w <-s_w + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
  title = element_text(color="black", size=30))
s_w

south_dry$ES <- south_dry$ES + 0.0001
south_dry$Trait <- factor(south_dry$Trait, level = c("SLA", "Leaf #", "Diameter", "Height", "Biomass","Germ","Flowering","Flower #", "Root"), labels = c("SLA", "Leaf #", "Diameter", "Height", "Biomass","Germ","Flowering","Flower #", "Root"))

s_d <- ggplot(south_dry, aes(Trait, y=ES, fill = Factor))+
  geom_bar(stat="identity", position = "dodge")+
  ggtitle("South Dry")+
  ylim(0, 0.19) +
  scale_color_manual(values= c("Year"="#FFA100", "TRTHIS"="skyblue3"))+
  scale_fill_manual(values= c("Year"="#FFA100", "TRTHIS"="skyblue3"))+
  theme_classic() +
  theme(legend.position = "none")
s_d <-s_d + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
  title = element_text(color="black", size=30))
s_d

library(cowplot)

gridExtra::grid.arrange(n_wet, n_dry, s_w, s_d)




north_wet %>%
  group_by(Factor) %>%
  summarise(
    count = n(),
    mean = mean(ES),
    sd = sd(ES),
  )
differences <- with(north_wet, ES[Factor == "Year"] - ES[Factor == "TRTHIS"])
shapiro.test(differences)
t.test(ES~Factor, data = north_wet, paired = TRUE)



north_dry %>%
  group_by(Factor) %>%
  summarise(
    count = n(),
    mean = mean(ES),
    sd = sd(ES),
  )
differences <- with(north_dry, ES[Factor == "Year"] - ES[Factor == "TRTHIS"])
shapiro.test(differences)
t.test(ES~Factor, data = north_dry, paired = TRUE)



south_wet %>%
  group_by(Year) %>%
  summarise(
    count = n(),
    mean = mean(ES),
    sd = sd(ES),
  )
differences <- with(south_wet, ES[Year == "Year"] - ES[Year == "TRTHIS"])
shapiro.test(differences)
wilcox.test(ES~Year, data = south_wet,paired=TRUE) #significantly different
t.test(ES~Year, data = south_wet, paired = TRUE)


south_dry %>%
  group_by(Factor) %>%
  summarise(
    count = n(),
    mean = mean(ES),
    sd = sd(ES),
  )
differences <- with(south_dry, ES[Factor == "Year"] - ES[Factor == "TRTHIS"])
shapiro.test(differences)
wilcox.test(ES~Factor, data = south_dry,paired=TRUE) #significantly different
