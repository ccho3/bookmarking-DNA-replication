# Generate plots from analysis1 results (related to Figure 1)

library(tidyverse)
source("custom_theme.R")

data<-read_csv("../results/Brd4_Cdc7_NC11-14.csv")

# Summarize data
summary <- data %>%
  group_by(Group, Frame, Channel) %>%
  summarize(
    mean_Intensity_Mean = mean(Intensity_Mean),
    SEM_Intensity_Mean = sd(Intensity_Mean) / sqrt(n()),
    mean_Intensity_Variance = mean(Intensity_Variance),
    SEM_Intensity_Variance = sd(Intensity_Variance) / sqrt(n())
    ) %>%
  mutate(
    Time = Frame*20,
    mean_Intensity_Variance = mean_Intensity_Variance / 100,
    SEM_Intensity_Variance = SEM_Intensity_Variance / 100
    ) %>%
  select(-Frame) %>%
  select(Time, everything())


# Rename Group
summary <- summary %>%
  mutate(Group = case_when(
    Group == "M10" ~ "NC11",
    Group == "M11" ~ "NC12",
    Group == "M12" ~ "NC13",
    Group == "M13" ~ "NC14"
  )) 

# Subset data
summary_Brd4 <- summary %>% filter(Channel=="Ch1")
summary_Cdc7 <- summary %>% filter(Channel=="Ch2")

# Set additional ggplot theme
theme_fig1<-theme(legend.position=c(0.5, 1.1),
                legend.direction = "horizontal",
                legend.text = element_text(size=6),
                legend.background = element_blank(),
                legend.title = element_blank(),
                legend.margin = margin(0, 0, 0, 0),
                legend.spacing = unit(0, "cm"),
                legend.key.size = unit(0.3, 'cm'))


# Generate plots
# svg, 200x165

# Brd4 mean
ggplot(summary_Brd4)+
  geom_line(aes(x=Time, y=mean_Intensity_Mean, group=Group, color=Group),linewidth=0.5)+
  geom_ribbon(aes(x=Time, 
                  ymin=mean_Intensity_Mean - SEM_Intensity_Mean, 
                  ymax=mean_Intensity_Mean + SEM_Intensity_Mean, 
                  group=Group, fill=Group), alpha=0.2)+
  scale_x_continuous(limits = c(0, 120), breaks=seq(0, 120, 20))+
  scale_y_continuous(limits = c(0, 170), breaks=seq(0, 200, 50))+
  theme_fig1

# Brd4 variance
ggplot(summary_Brd4)+
  geom_line(aes(x=Time, y=mean_Intensity_Variance, group=Group, color=Group),linewidth=0.5)+
  geom_ribbon(aes(x=Time, 
                  ymin=mean_Intensity_Variance - SEM_Intensity_Variance, 
                  ymax=mean_Intensity_Variance + SEM_Intensity_Variance, 
                  group=Group, fill=Group), alpha=0.2)+
  scale_x_continuous(limits = c(0, 120), breaks=seq(0, 120, 20))+
  scale_y_continuous(limits = c(0, 330), breaks=seq(0, 300, 50))+
  theme_fig1

# Cdc7 mean
ggplot(summary_Cdc7)+
  geom_line(aes(x=Time, y=mean_Intensity_Mean, group=Group, color=Group),linewidth=0.5)+
  geom_ribbon(aes(x=Time, 
                  ymin=mean_Intensity_Mean - SEM_Intensity_Mean, 
                  ymax=mean_Intensity_Mean + SEM_Intensity_Mean, 
                  group=Group, fill=Group), alpha=0.2)+
  scale_x_continuous(limits = c(0, 120), breaks=seq(0, 120, 20))+
  scale_y_continuous(limits = c(0, 550), breaks=seq(0, 500, 100))+
  theme_fig1

# Cdc7 variance
ggplot(summary_Cdc7)+
  geom_line(aes(x=Time, y=mean_Intensity_Variance, group=Group, color=Group),linewidth=0.5)+
  geom_ribbon(aes(x=Time, 
                  ymin=mean_Intensity_Variance - SEM_Intensity_Variance, 
                  ymax=mean_Intensity_Variance + SEM_Intensity_Variance, 
                  group=Group, fill=Group), alpha=0.2)+
  scale_x_continuous(limits = c(0, 120), breaks=seq(0, 120, 20))+
  scale_y_continuous(limits = c(0, 850), breaks=seq(0, 800, 200))+
  theme_fig1

