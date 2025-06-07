# Generate plots from analysis1 results (related to Figure 1)

library(tidyverse)
source("custom_theme.R")

data<-read_csv("../results/Brd4_Cdc7_NC11-14.csv")

# Get peak mean or variance for each movie

data_peak_mean <- data %>%
  group_by(Group, Replicate, Channel) %>%
  slice_max(order_by = Intensity_Mean) %>%
  arrange(Channel) %>%
  select(-Intensity_Variance) %>%
  mutate(Group = case_when(
    Group == "M10" ~ "NC11",
    Group == "M11" ~ "NC12",
    Group == "M12" ~ "NC13",
    Group == "M13" ~ "NC14"
  ))

data_peak_variance <- data %>%
  group_by(Group, Replicate, Channel) %>%
  slice_max(order_by = Intensity_Variance) %>%
  arrange(Channel) %>%
  select(-Intensity_Mean) %>%
  mutate(Group = case_when(
    Group == "M10" ~ "NC11",
    Group == "M11" ~ "NC12",
    Group == "M12" ~ "NC13",
    Group == "M13" ~ "NC14"
  )) %>%
  mutate(
    Intensity_Variance = Intensity_Variance / 100
    )


summary_peak_mean <- data_peak_mean %>%
  group_by(Group, Channel) %>%
  summarize(
    mean_Peak_Mean = mean(Intensity_Mean),
    SEM_Peak_Mean = sd(Intensity_Mean) / sqrt(n())
  ) %>%
  arrange(Channel)
 
summary_peak_variance <- data_peak_variance %>%
  group_by(Group, Channel) %>%
  summarize(
    mean_Peak_Variance = mean(Intensity_Variance),
    SEM_Peak_Variance = sd(Intensity_Variance) / sqrt(n())
  ) %>%
  arrange(Channel) 



# Set additional ggplot theme
theme_fig1<-theme(legend.position=c(0.5, 1.1),
                  legend.direction = "horizontal",
                  legend.text = element_text(size=6),
                  legend.background = element_blank(),
                  legend.title = element_blank(),
                  legend.margin = margin(0, 0, 0, 0),
                  legend.spacing = unit(0, "cm"),
                  legend.key.size = unit(0.3, 'cm'))



# Brd4 mean, 200x200
ggplot(summary_peak_mean %>% filter(Channel == "Ch1"), 
       aes(x = Group, y = mean_Peak_Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(
    data = data_peak_mean %>% filter(Channel == "Ch1"),
    aes(x = Group, y = Intensity_Mean), size = 0.5, color = "black", alpha = 0.6,
    show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean_Peak_Mean - SEM_Peak_Mean,
        ymax = mean_Peak_Mean + SEM_Peak_Mean),
    width = 0.2,
    position = position_dodge(width = 0.8),
    linewidth = 0.3
  ) +
  theme(legend.position = "none")


# Brd4 variance, 200x200
ggplot(summary_peak_variance %>% filter(Channel == "Ch1"), 
       aes(x = Group, y = mean_Peak_Variance, fill = Group)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(
    data = data_peak_variance %>% filter(Channel == "Ch1"),
    aes(x = Group, y = Intensity_Variance), size = 0.5, color = "black", alpha = 0.6,
    show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean_Peak_Variance - SEM_Peak_Variance,
        ymax = mean_Peak_Variance + SEM_Peak_Variance),
    width = 0.2,
    position = position_dodge(width = 0.8),
    linewidth = 0.3
  ) +
  theme(legend.position = "none")



# Cdc7 mean, 200x200
ggplot(summary_peak_mean %>% filter(Channel == "Ch2"), 
       aes(x = Group, y = mean_Peak_Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(
    data = data_peak_mean %>% filter(Channel == "Ch2"),
    aes(x = Group, y = Intensity_Mean), size = 0.5, color = "black", alpha = 0.6,
    show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean_Peak_Mean - SEM_Peak_Mean,
        ymax = mean_Peak_Mean + SEM_Peak_Mean),
    width = 0.2,
    position = position_dodge(width = 0.8),
    linewidth = 0.3
  ) +
  theme(legend.position = "none")


# Cdc7 variance, 200x200
ggplot(summary_peak_variance %>% filter(Channel == "Ch2"), 
       aes(x = Group, y = mean_Peak_Variance, fill = Group)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(
    data = data_peak_variance %>% filter(Channel == "Ch2"),
    aes(x = Group, y = Intensity_Variance), size = 0.5, color = "black", alpha = 0.6,
    show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = mean_Peak_Variance - SEM_Peak_Variance,
        ymax = mean_Peak_Variance + SEM_Peak_Variance),
    width = 0.2,
    position = position_dodge(width = 0.8),
    linewidth = 0.3
  ) +
  theme(legend.position = "none")















