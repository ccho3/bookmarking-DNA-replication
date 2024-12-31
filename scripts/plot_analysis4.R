# Generate plot from analysis4 results (related to Figure 4)

# Load package and data
library(tidyverse)

data<-read_csv("../results/Rpb3_MCP_summary.csv")

# Calculate percentages of nuclei
# Convert frame to time in minutes
data <- data %>%
  mutate(
    cumulative_MCP_ON = cumulative_MCP_ON*100 / total_nuclei,
    cumulative_MCP_off = cumulative_MCP_off*100 / total_nuclei,
    Time = 2 + frame*0.5,
    Treatment = case_when(
      treatment == "water" ~ "Control",
      treatment == "XL413" ~ "Cdc7i"),
    Treatment = factor(Treatment, levels=c("Control", "Cdc7i"))
    )

# Add data points for Time = 0-1.5 min
# Images were visually inspected to confirm the absence of MCP foci
time_points <- tibble(Time = seq(0, 1.5, 0.5))

blank_data <- data %>%
  select(Treatment, cycle, replicate, total_nuclei) %>%
  distinct() %>%
  full_join(time_points, by = character()) %>%
  mutate(
    cumulative_MCP_ON = 0,
    cumulative_MCP_off = 0
  )

# Bind the blank data rows to the original data
data2 <- bind_rows(data, blank_data) %>%
  arrange(Treatment, cycle, replicate, Time) %>%
  select(Treatment, cycle, replicate, Time, everything())

# Rearrange data for plotting
data_long <- data2 %>%
  pivot_longer(cols = c(cumulative_MCP_ON, cumulative_MCP_off), 
               names_to = "MCP foci", 
               names_prefix = "cumulative_MCP_",
               values_to = "Percentage") %>%
  mutate(
    `MCP foci` = case_when(
      `MCP foci` == "ON" ~ "Appeared",
      `MCP foci` == "off" ~ "Disappeared"),
    `MCP foci` = factor(`MCP foci`, levels=c("Appeared", "Disappeared"))
  ) %>%
  select(-c(instantaneous_MCP_ON, total_nuclei, treatment, frame)) %>%
  arrange(`MCP foci`)

# Calculate mean and SEM
data_long_summary <- data_long %>%
  group_by(Treatment, cycle, Time, `MCP foci`) %>%
  summarize(
    mean_Percentage = mean(Percentage),
    SEM_Percentage = sd(Percentage) / sqrt(n())
    ) %>%
  arrange(`MCP foci`)

# Adjust Time to be continuous across cycles
data_long_summary <- data_long_summary %>%
  mutate(
    Time_cumulative = case_when(
      cycle=="nc12" ~ Time,
      cycle=="nc13" ~ Time +12,
      cycle=="nc14" ~ Time +31.5)
    
    )

# Create custome breaks and labels for x-axis
x_breaks <- c(seq(0, 7, 2), seq(0, 15, 2)+12, seq(0, 16, 2)+31.5)
x_labels <- c(seq(0, 7, 2), seq(0, 15, 2), seq(0, 16, 2))

##### Set ggplot theme
my_theme<-theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.background = element_blank(), 
                panel.border = element_rect(fill=NA, linewidth=0.5),
                text = element_text(size=7),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text = element_text(size=7, colour="black"),
                axis.line = element_blank(),
                axis.ticks = element_line(linewidth =0.5),
                plot.margin = unit(c(0,0,0,0),"cm"),
                legend.position = "right",
                legend.margin = margin(0, 0, 0, 0),
                legend.spacing = unit(0.1, "cm"),
                legend.spacing.y = unit(0.4, "cm"),
                legend.key.size = unit(0.5, "cm"),
                legend.text = element_text(size=7)
                )

theme_set(my_theme)

# Plot
# 780x190
ggplot(data_long_summary)+
  geom_line(aes(x=Time_cumulative, y=mean_Percentage, group=interaction(Treatment, `MCP foci`, cycle),
                color=Treatment, linetype=`MCP foci`),linewidth=0.5)+
  geom_ribbon(aes(x = Time_cumulative, ymin = mean_Percentage - SEM_Percentage, 
                  ymax = mean_Percentage + SEM_Percentage, 
                  fill = Treatment, group = interaction(Treatment, `MCP foci`, cycle)), 
              alpha = 0.2)+
  annotate("rect", xmin=7, xmax=12, ymin=-Inf, ymax=Inf, fill="lightgrey", alpha=0.5)+
  annotate("rect", xmin=26.5, xmax=31.5, ymin=-Inf, ymax=Inf, fill ="lightgrey", alpha=0.5)+
  scale_color_manual(values=c("#619CFF","#F8766D"))+
  scale_fill_manual(values=c("#619CFF","#F8766D"))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_x_continuous(limits = c(-1, 48), labels=x_labels, breaks=x_breaks, expand=c(0, 0))+
  scale_y_continuous(limits = c(-2, 100), breaks=seq(0, 100, 25))+
  guides(color = guide_legend(order = 1),  # Treatment (color) legend first
         fill = guide_legend(order = 1),   # Combine color and fill legends
         linetype = guide_legend(order = 2))




