#R 3.5.1
library(ggplot2)
library(gridExtra)
#Slightly edited theme for heatmap
ggplot_theme <-theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=16),
        axis.title.x   = element_text(size=16, face="bold"),
        axis.title.y   = element_text(size=16, face="bold"),
        strip.text.x = element_text(size=12, face="bold"),
        strip.background = element_blank(),
        legend.text = element_text(size = 13),
        #legend.title =  element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank()
        #axis.line.x = element_line(colour = "black"), 
        #axis.line.y = element_line(colour = "black")
  )
theme_set(ggplot_theme)

#Load data
sigFishers_enriched_subfamilies_enrichBases_forR <- read.delim("sigFishers_enriched_subfamilies_enrichBases_forR")

#Plot
ggplot(sigFishers_enriched_subfamilies_enrichBases_forR, aes(Subfamily, cCRE, fill=FoldEnrich), group=Class) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='gray') +
  labs(x="", y="", fill="log2\nenrichment") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid(~Class, scales='free_x', space="free_x") #Save as 8x30 landscape

