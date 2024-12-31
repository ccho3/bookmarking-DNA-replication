# Set ggplot theme

my_theme<-theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(fill=NA, linewidth=0.5),
                plot.background = element_blank(),
                plot.margin=unit(c(0.75,0,0,0),"cm"),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text = element_text(size=6, colour="black"),
                axis.line = element_blank(),
                axis.ticks = element_line(linewidth =0.5),
                text = element_text(size=6)
                )

theme_set(my_theme)
