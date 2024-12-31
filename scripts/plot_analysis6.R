# Generate plot from analysis6 results (related to Figure S3)

library(tidyverse)
source("custom_theme.R")

data<-read_csv("../results/PCNA_Rpb3_max_Cdc7i.csv")

# Summarize data
data2 <- data %>%
  group_by(Group, Frame, Channel) %>%
  summarize(
    mean = mean(Intensity_Max),
    SEM = sd(Intensity_Max)/sqrt(n())
    ) %>%
  mutate(
    Time = Frame*0.5,
    mean = mean / 10,
    SEM = SEM / 10
    ) %>%
  select(-Frame) %>%
  select(Time, everything()) %>%
  mutate(
    Group=case_when(
      Group == "XL413" ~ "Cdc7i",
      Group == "water" ~ "Water")
  ) %>% 
  filter(Time <= 6)


#####

# Set additional ggplot theme
theme1<-theme(legend.position=c(0.5, 1.1),
                  legend.direction = "horizontal",
                  legend.text = element_text(size=6),
                  legend.background = element_blank(),
                  legend.margin = margin(0, 0, 0, 0),
                  legend.spacing = unit(0, "cm"),
                  legend.title = element_blank(),
                  legend.key.size = unit(0.3, 'cm'))


# 200x165
# Max intensity
ggplot(data2%>%filter(Channel=="Ch2"))+
  geom_line(aes(x=Time, y=mean, group=factor(Group, levels=c("Water", "Cdc7i")),
                color=factor(Group, levels=c("Water", "Cdc7i"))),linewidth=0.5)+
  geom_ribbon(aes(x=Time, 
                  ymin=mean - SEM, 
                  ymax=mean + SEM, 
                  group=factor(Group, levels=c("Water", "Cdc7i")), 
                  fill=factor(Group, levels=c("Water", "Cdc7i"))), alpha=0.2)+
  scale_y_continuous(limits = c(0, 1000), breaks=seq(0, 1000, 200))+
  scale_x_continuous(limits = c(0, 6), breaks=seq(0, 6, 1))+
  scale_color_manual(values=c("#619CFF","#F8766D"))+
  scale_fill_manual(values=c("#619CFF","#F8766D"))+
  theme1

