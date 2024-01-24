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

#Subfamily enrichments first
#Load data
sigFishers_subfamilies_enrichBases_forR <- read.delim("sigFishers_subfamilies_enrichBases_forR")

#Plot (Don't use this, doesn't work well when including depleted subfamilies)
ggplot(sigFishers_subfamilies_enrichBases_forR, aes(Subfamily, cCRE, fill=FoldEnrich), group=Class) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='gray') +
  labs(x="", y="", fill="log2\nenrichment") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_wrap(~Class, scales='free_x') #Save as 20x40 landscape

#Class level enrichments
#Load data
classLevel_cCRE_enrichments_all_combined_forR <- read.delim("classLevel_cCRE_enrichments_all_combined_forR")
classLevel_cCRE_enrichments_all_combined_forR$cCRE <- 
  factor(classLevel_cCRE_enrichments_all_combined_forR$cCRE, 
         levels=c("dELS", "pELS", "PLS", "DNase-H3K4me3", "CTCF-only"))
classLevel_cCRE_enrichments_all_combined_forR <- subset(classLevel_cCRE_enrichments_all_combined_forR,
                                                        Class != "SVA" & Class != "Other")
#Plot
ggplot(classLevel_cCRE_enrichments_all_combined_forR, aes(cCRE, Class, fill=FoldEnrich)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='gray') +
  labs(x="", y="", fill="log2\nenrichment", title="TE class overlap with ENCODE cCREs") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#Try boxplots for distribution of subfamilies within classes
#Load data
subfamily_cCRE_enrichments_all_combined_forR <- read.delim("subfamily_cCRE_enrichments_all_combined_forR")
subfamily_cCRE_enrichments_all_combined_forR$cCRE <- 
  factor(subfamily_cCRE_enrichments_all_combined_forR$cCRE, 
         levels=c("dELS", "pELS", "PLS", "DNase-H3K4me3", "CTCF-only"))

##Use this only if including all 0 overlap subfamilies (NA as enrichment value)
#Convert "NA" enrichment values to -10 (pseudocount)
#library(tidyr)
#subfamily_cCRE_enrichments_all_combined_forR$FoldEnrich <- 
#  replace_na(subfamily_cCRE_enrichments_all_combined_forR$FoldEnrich, -10)

#Plot with SVA, horizontal bars
ggplot(subfamily_cCRE_enrichments_all_combined_forR, aes(FoldEnrich, cCRE, fill=Class)) +
  geom_boxplot() +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(x="log2 fold enrichment", y="", title="TE subfamily ENCODE cCRE enrichment distribution by class") +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank()
        ,panel.grid.major.x = element_line(colour="gray")
  )
#Plot with SVA, vertical bars
ggplot(subfamily_cCRE_enrichments_all_combined_forR, aes(cCRE, FoldEnrich, fill=Class)) +
  geom_boxplot() +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(y="log2 fold enrichment", x="", title="TE subfamily ENCODE cCRE enrichment distribution by class") +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_line(colour="gray"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )
#Plot without SVA, horizontal bars
subfamily_cCRE_enrichments_all_combined_forR_noSVA <- subset(subfamily_cCRE_enrichments_all_combined_forR, Class != "SVA")
ggplot(subfamily_cCRE_enrichments_all_combined_forR_noSVA, aes(FoldEnrich, cCRE, fill=Class)) +
  geom_boxplot() +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(x="log2 fold enrichment", y="", title="TE subfamily ENCODE cCRE enrichment distribution by class") +
  scale_fill_manual("TE class", values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")) +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank()
        ,panel.grid.major.x = element_line(colour="gray")
  )
#Plot without SVA, vertical bars
ggplot(subfamily_cCRE_enrichments_all_combined_forR_noSVA, aes(cCRE, FoldEnrich, fill=Class)) +
  geom_boxplot() +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(y="log2 fold enrichment", x="", title="TE subfamily ENCODE cCRE enrichment distribution by class") +
  scale_fill_manual("TE class", values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")) +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_line(colour="gray"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

#Subfamily-cCRE enrichments for 25 full classified cell/tissue types
#Use minimum overlap threshold of 50 for subfamilies to include
#Load data
enrichments_TEsubfamilies_fullClassCCREs_dELS <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_dELS")
enrichments_TEsubfamilies_fullClassCCREs_pELS <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_pELS")
enrichments_TEsubfamilies_fullClassCCREs_PLS <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_PLS")
enrichments_TEsubfamilies_fullClassCCREs_DNase.H3K4me3 <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_DNase-H3K4me3")
enrichments_TEsubfamilies_fullClassCCREs_DNase.only <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_DNase-only")
enrichments_TEsubfamilies_fullClassCCREs_CTCF.only <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_CTCF-only")
enrichments_TEsubfamilies_fullClassCCREs_CTCF.bound <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_CTCF-bound")
#Plot
ggplot(enrichments_TEsubfamilies_fullClassCCREs_dELS, aes(Subfamily, Celltype, fill=Enrichment)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='grey') +
  #scale_fill_gradient("viridis") +
  labs(x="", y="", fill="log2\nenrichment", title="TE class overlap with ENCODE cCREs - dELS") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggplot(enrichments_TEsubfamilies_fullClassCCREs_pELS, aes(Subfamily, Celltype, fill=Enrichment)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='grey') +
  #scale_fill_gradient("viridis") +
  labs(x="", y="", fill="log2\nenrichment", title="TE class overlap with ENCODE cCREs - pELS") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggplot(enrichments_TEsubfamilies_fullClassCCREs_PLS, aes(Subfamily, Celltype, fill=Enrichment)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='grey') +
  #scale_fill_gradient("viridis") +
  labs(x="", y="", fill="log2\nenrichment", title="TE class overlap with ENCODE cCREs - PLS") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggplot(enrichments_TEsubfamilies_fullClassCCREs_DNase.H3K4me3, aes(Subfamily, Celltype, fill=Enrichment)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='grey') +
  #scale_fill_gradient("viridis") +
  labs(x="", y="", fill="log2\nenrichment", title="TE class overlap with ENCODE cCREs - DNase-H3K4me3") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggplot(enrichments_TEsubfamilies_fullClassCCREs_DNase.only, aes(Subfamily, Celltype, fill=Enrichment)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='grey') +
  #scale_fill_gradient("viridis") +
  labs(x="", y="", fill="log2\nenrichment", title="TE class overlap with ENCODE cCREs - DNase only") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggplot(enrichments_TEsubfamilies_fullClassCCREs_CTCF.only, aes(Subfamily, Celltype, fill=Enrichment)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='grey') +
  #scale_fill_gradient("viridis") +
  labs(x="", y="", fill="log2\nenrichment", title="TE class overlap with ENCODE cCREs - CTCF only") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggplot(enrichments_TEsubfamilies_fullClassCCREs_CTCF.bound, aes(Subfamily, Celltype, fill=Enrichment)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", na.value='grey') +
  #scale_fill_gradient("viridis") +
  labs(x="", y="", fill="log2\nenrichment", title="TE class overlap with ENCODE cCREs - all CTCF bound") +
  theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#Subfamily-cCRE enrichments for 25 full classified cell/tissue types, separate by TE class
#Use minimum overlap threshold of 10 for subfamilies to include
#Load data
enrichments_TEsubfamilies_fullClassCCREs_wClass_dELS <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass_dELS")
enrichments_TEsubfamilies_fullClassCCREs_wClass_pELS <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass_pELS")
enrichments_TEsubfamilies_fullClassCCREs_wClass_PLS <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass_PLS")
enrichments_TEsubfamilies_fullClassCCREs_wClass_DNase.H3K4me3 <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass_DNase-H3K4me3")
enrichments_TEsubfamilies_fullClassCCREs_wClass_DNase.only <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass_DNase-only")
enrichments_TEsubfamilies_fullClassCCREs_wClass_CTCF.only <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass_CTCF-only")
enrichments_TEsubfamilies_fullClassCCREs_wClass_CTCF.bound <- read.delim("forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass_CTCF-bound")
#Function for plotting TE classes separately and then combining
plot_class_enrich_separate <- function(input_df) {
  scale_max <- max(input_df$Enrichment, na.rm = TRUE)
  scale_min <- min(input_df$Enrichment, na.rm = TRUE)
  p1 <- ggplot(input_df[input_df$Class=='DNA',], aes(Subfamily, Celltype, fill=Enrichment)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max)) +
    labs(x="", y="", fill="log2\nenrichment", title="DNA") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p2 <- ggplot(input_df[input_df$Class=='LINE',], aes(Subfamily, Celltype, fill=Enrichment)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max)) +
    labs(x="", y="", fill="log2\nenrichment", title="LINE") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p3 <- ggplot(input_df[input_df$Class=='SINE',], aes(Subfamily, Celltype, fill=Enrichment)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max)) +
    labs(x="", y="", fill="log2\nenrichment", title="SINE") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p4 <- ggplot(input_df[input_df$Class=='LTR',], aes(Subfamily, Celltype, fill=Enrichment)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max)) +
    labs(x="", y="", fill="log2\nenrichment", title="LTR") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p5 <- ggplot(input_df[input_df$Class=='SVA',], aes(Subfamily, Celltype, fill=Enrichment)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max)) +
    labs(x="", y="", fill="log2\nenrichment", title="SVA") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  lay <- rbind(c(1, 1, 1, 1, NA, NA), c(2, 2, 3, 5, NA, NA), c(4, 4, 4, 4, 4, 4))
  grid.arrange(p1, p2, p3, p4, p5, ncol=1, layout_matrix=lay)
}
#Plot (save as 30x100?)
plot_class_enrich_separate(enrichments_TEsubfamilies_fullClassCCREs_wClass_dELS)
plot_class_enrich_separate(enrichments_TEsubfamilies_fullClassCCREs_wClass_pELS)
plot_class_enrich_separate(enrichments_TEsubfamilies_fullClassCCREs_wClass_PLS)
plot_class_enrich_separate(enrichments_TEsubfamilies_fullClassCCREs_wClass_DNase.H3K4me3)
plot_class_enrich_separate(enrichments_TEsubfamilies_fullClassCCREs_wClass_DNase.only)
plot_class_enrich_separate(enrichments_TEsubfamilies_fullClassCCREs_wClass_CTCF.only)
plot_class_enrich_separate(enrichments_TEsubfamilies_fullClassCCREs_wClass_CTCF.bound)

