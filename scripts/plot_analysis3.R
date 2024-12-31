# Generate plots from analysis3 results (related to Figure 3)

library(tidyverse)
source("custom_theme.R")

data<-read_csv("../results/Brd4_Cdc7_max_Cdc7i.csv")

# Transform data
data2 <- data %>%
  mutate(max_Brd4 = Brd4_max/3,
         max_Cdc7 = Cdc7_max/4,
         Treatment = case_when(
           Treatment == "water" ~ "Water",
           Treatment == "XL413" ~ "Cdc7i"),
         Time = case_when(
           Time == "0sec" ~ "0 sec",
           Time == "80sec" ~ "80 sec")
         ) %>%
  pivot_longer(cols = starts_with("max_"),
               names_to = "Channel",
               names_prefix = "max_",
               values_to = "Max"
  ) %>%
  select(Treatment, Replicate, Time, Channel, Max)

# Subset data
data_0sec <- data2 %>% filter(Time=="0 sec") 
data_80sec <- data2 %>% filter(Time=="80 sec") 

# Set ggplot theme
theme_fig3<-theme(text = element_text(size=8),
                  plot.margin=unit(c(0.1, 0.7, 0.7, 0.7),"cm"),
                  axis.text = element_text(size=7, colour="black"),
                  legend.position="bottom",
                  legend.text = element_text(size=7),
                  legend.margin = margin(0, 0, 0, 0),
                  legend.background = element_blank(),
                  legend.title = element_blank(),
                  legend.key.size = unit(0.4, 'cm'))

# Generate plots
# svg, 240x220

ggplot(data_0sec, aes(x=Channel, y=Max, 
                  fill=factor(Treatment, levels=c("Water", "Cdc7i"))))+
  geom_boxplot(size=0.3, outlier.size = 0.5)+
  scale_y_continuous(limits = c(0, 1100), breaks=seq(0, 1000, 200))+
  scale_fill_manual(values=c("#619CFF","#F8766D"))+
  geom_segment(aes(x = 0.75, y = 1000, xend = 1.25, yend = 1000), linewidth=0.3)+
  geom_segment(aes(x = 1.75, y = 1000, xend = 2.25, yend = 1000), linewidth=0.3)+
  annotate("text", x = 1, y = 1020, label= "", size=4)+ 
  annotate("text", x=2, y = 1020, label= "", size=4)+
  theme_fig3

ggplot(data_80sec, aes(x=Channel, y=Max, 
                      fill=factor(Treatment, levels=c("Water", "Cdc7i"))))+
  geom_boxplot(size=0.3, outlier.size = 0.5)+
  scale_y_continuous(limits = c(0, 1100), breaks=seq(0, 1000, 200))+
  scale_fill_manual(values=c("#619CFF","#F8766D"))+
  geom_segment(aes(x = 0.75, y = 1000, xend = 1.25, yend = 1000), linewidth=0.3)+
  geom_segment(aes(x = 1.75, y = 1000, xend = 2.25, yend = 1000), linewidth=0.3)+
  annotate("text", x = 1, y = 1020, label= "", size=4)+
  annotate("text", x=2, y = 1020, label= "", size=4)+
  theme_fig3


# Mann-Whitney U test
p_values <- c(
  wilcox.test(Max ~ Treatment, data_0sec %>% filter(Channel=="Brd4"))$p.value,
  wilcox.test(Max ~ Treatment, data_0sec %>% filter(Channel=="Cdc7"))$p.value,
  wilcox.test(Max ~ Treatment, data_80sec %>% filter(Channel=="Brd4"))$p.value,
  wilcox.test(Max ~ Treatment, data_80sec %>% filter(Channel=="Cdc7"))$p.value
)

adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
adjusted_p_values


#
sample_size <- data2 %>% 
  count(Treatment, Time, Channel) %>%
  arrange(Time, Channel, Treatment)

sample_size

