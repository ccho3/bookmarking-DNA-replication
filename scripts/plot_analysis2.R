# Generate plots from analysis2 results (related to Figure 2 and Figure S2)

library(tidyverse)
library(gridExtra)
source("custom_theme.R")

data<-read_csv("../results/Brd4_Cdc7_treatments.csv")

# Prepare data
data2 <- data %>%
  mutate(Cdc7_variance = Cdc7_variance/1000,
         Experiments = case_when(
           Group == "DMSO" ~ "Injection",
           Group == "JQ1" ~ "Injection",
           Group == "w" ~ "RNAi",
           Group == "nej" ~ "RNAi"),
         Brd4_mean = ifelse(Brd4_mean < 0, 0, Brd4_mean)
         )
  


# Set additional ggplot theme
theme_plot2<-theme(panel.border = element_rect(fill=NA, linewidth=0.75),
                   text = element_text(size=8),
                   axis.text = element_text(size=8, colour="black"),
                   axis.ticks = element_line(linewidth =0.55),
                   plot.margin=unit(c(0.7, 0.7, 0.7, 0.7),"cm"),
                   legend.position=c(0.15, 0.9),
                   legend.text = element_text(size=8),
                   legend.background = element_blank(),
                   legend.title = element_blank(),
                   legend.key.size = unit(0.3, 'cm')
                   )


# Plots
# Brd4 intensity
height1<-280
line_pos1<-266
a<-ggplot(data2, aes(x=factor(Group, levels = c("w", "nej", "DMSO", "JQ1")), 
                          y=Brd4_mean))+
  geom_boxplot(size=0.4, outlier.size = 0.6)+
  scale_y_continuous(limits = c(-14, height1), breaks=seq(0, height1, 100))+
  geom_segment(x = 1, y = line_pos1, xend = 2, yend = line_pos1)+
  geom_segment(x = 3, y = line_pos1, xend = 4, yend = line_pos1)+
  annotate("text", x=1.5, y=line_pos1+5, label= "***", size=5)+
  annotate("text", x=3.5, y=line_pos1+5, label= "***", size=5)+
  theme_plot2
  
# Cdc7 intensity
height2<-530
line_pos2<-515
b<-ggplot(data2, aes(x=factor(Group, levels = c("w", "nej", "DMSO", "JQ1")), 
                  y=Cdc7_mean))+
  geom_boxplot(size=0.4, outlier.size = 0.6)+
  scale_y_continuous(limits = c(-27, height2), breaks=seq(0, height2, 100))+
  geom_segment(x = 1, y = line_pos2, xend = 2, yend = line_pos2)+
  geom_segment(x = 3, y = line_pos2, xend = 4, yend = line_pos2)+
  annotate("text", x=1.5, y=line_pos2+5, label= "***", size=5)+
  annotate("text", x=3.5, y=line_pos2+5, label= "***", size=5)+
  theme_plot2

# Cdc7 variance
height3<-195
line_pos3<-185
c<-ggplot(data2, aes(x=factor(Group, levels = c("w", "nej", "DMSO", "JQ1")), 
                 y=Cdc7_variance))+
  geom_boxplot(size=0.4, outlier.size = 0.6)+
  scale_y_continuous(limits = c(-10, height3), breaks=seq(0, height3, 100))+
  geom_segment(x = 1, y = line_pos3, xend = 2, yend = line_pos3)+
  geom_segment(x = 3, y = line_pos3, xend = 4, yend = line_pos3)+
  annotate("text", x=1.5, y=line_pos3+5, label= "***", size=5)+
  annotate("text", x=3.5, y=line_pos3+5, label= "***", size=5)+
  theme_plot2

grid.arrange(a, b, c, ncol=3)

#
p_values <- c(
  wilcox.test(Brd4_mean ~ Group, data2 %>% filter(Experiments == "RNAi"))$p.value,
  wilcox.test(Brd4_mean ~ Group, data2 %>% filter(Experiments == "Injection"))$p.value,
  wilcox.test(Cdc7_mean ~ Group, data2 %>% filter(Experiments == "RNAi"))$p.value,
  wilcox.test(Cdc7_mean ~ Group, data2 %>% filter(Experiments == "Injection"))$p.value,
  wilcox.test(Cdc7_variance ~ Group, data2 %>% filter(Experiments == "RNAi"))$p.value,
  wilcox.test(Cdc7_variance ~ Group, data2 %>% filter(Experiments == "Injection"))$p.value
)

adjusted_p_values <- p.adjust(p_values, method = "bonferroni")


#
sample_size <- data2 %>% count(Group)
sample_size
