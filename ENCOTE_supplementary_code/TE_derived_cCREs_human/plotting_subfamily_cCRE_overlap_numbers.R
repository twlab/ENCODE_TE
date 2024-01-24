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

#Distribution of cCRE overlaps (number) in TE subfamilies separated by class
#Load data
subfamily_cCRE_overlap_numbers_all_combined_forR_count <- read.delim("subfamily_cCRE_overlap_numbers_all_combined_forR_count")
subfamily_cCRE_overlap_numbers_all_combined_forR_count$cCRE <- 
  factor(subfamily_cCRE_overlap_numbers_all_combined_forR_count$cCRE, 
         levels=c("dELS", "pELS", "PLS", "DNase-H3K4me3", "CTCF-only"))
#Add pseudocount of 0.1 for subfamilies with 0 overlap
subfamily_cCRE_overlap_numbers_all_combined_forR_count$Number <- subfamily_cCRE_overlap_numbers_all_combined_forR_count$Number + 0.1
#Plot with SVA, horizontal bars
ggplot(subfamily_cCRE_overlap_numbers_all_combined_forR_count, aes(Number, cCRE, fill=Class)) +
  geom_boxplot() +
  scale_x_continuous(trans='log10') +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(x="Number of elements", y="", title="TE subfamily ENCODE cCRE overlap number distribution by class") +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank()
        ,panel.grid.major.x = element_line(colour="gray")
  )
#Plot with SVA, vertical bars
ggplot(subfamily_cCRE_overlap_numbers_all_combined_forR_count, aes(cCRE, Number, fill=Class)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(y="Number of elements", x="", title="TE subfamily ENCODE cCRE overlap number distribution by class") +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_line(colour="gray"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )
#Plot without SVA
subfamily_cCRE_overlap_numbers_all_combined_forR_count_noSVA <- subset(subfamily_cCRE_overlap_numbers_all_combined_forR_count, Class != "SVA")
ggplot(subfamily_cCRE_overlap_numbers_all_combined_forR_count_noSVA, aes(Number, cCRE, fill=Class)) +
  geom_boxplot() +
  scale_x_continuous(trans='log10') +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(x="Number of elements", y="", title="TE subfamily ENCODE cCRE overlap number distribution by class") +
  scale_fill_manual("TE class", values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")) +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank()
        ,panel.grid.major.x = element_line(colour="gray")
  )
#Plot without SVA, vertical bars
ggplot(subfamily_cCRE_overlap_numbers_all_combined_forR_count_noSVA, aes(cCRE, Number, fill=Class)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(y="Number of elements", x="", title="TE subfamily ENCODE cCRE overlap number distribution by class") +
  scale_fill_manual("TE class", values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")) +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank(),
        panel.grid.major.y = element_line(colour="gray"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

#Distribution of cCRE overlaps (base pairs) in TE subfamilies separated by class
#Load data
subfamily_cCRE_overlap_numbers_all_combined_forR_bases <- read.delim("subfamily_cCRE_overlap_numbers_all_combined_forR_bases")
#Add pseudocount of 0.1 for subfamilies with 0 overlap
subfamily_cCRE_overlap_numbers_all_combined_forR_bases$Number <- subfamily_cCRE_overlap_numbers_all_combined_forR_bases$Number + 0.1
#Plot with SVA
ggplot(subfamily_cCRE_overlap_numbers_all_combined_forR_bases, aes(Number, cCRE, fill=Class)) +
  geom_boxplot() +
  scale_x_continuous(trans='log10') +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(x="Bp of overlap", y="", title="TE subfamily ENCODE cCRE overlap distribution by class") +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank()
        ,panel.grid.major.x = element_line(colour="gray")
  )
#Plot without SVA
subfamily_cCRE_overlap_numbers_all_combined_forR_bases_noSVA <- subset(subfamily_cCRE_overlap_numbers_all_combined_forR_bases, Class != "SVA")
ggplot(subfamily_cCRE_overlap_numbers_all_combined_forR_bases_noSVA, aes(Number, cCRE, fill=Class)) +
  geom_boxplot() +
  scale_x_continuous(trans='log10') +
  #geom_vline(xintercept = 0, colour="gray") +
  labs(x="Bp of overlap", y="", title="TE subfamily ENCODE cCRE overlap distribution by class") +
  scale_fill_manual("TE class", values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")) +
  theme(axis.line.x = element_line(colour='black'), axis.line.y = element_line(colour='black'),
        axis.ticks = element_blank()
        ,panel.grid.major.x = element_line(colour="gray")
  )

#Subfamily-cCRE enrichments for 25 full classified cell/tissue types, separate by TE class
#Use minimum overlap threshold of 10 for subfamilies to include
#Load data
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_dELS <- read.delim("forR_plotting/overlapCounts_TEsubfamilies_fullClassCCREs_wClass_dELS")
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_pELS <- read.delim("forR_plotting/overlapCounts_TEsubfamilies_fullClassCCREs_wClass_pELS")
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_PLS <- read.delim("forR_plotting/overlapCounts_TEsubfamilies_fullClassCCREs_wClass_PLS")
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase.H3K4me3 <- read.delim("forR_plotting/overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase-H3K4me3")
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase.only <- read.delim("forR_plotting/overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase-only")
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF.only <- read.delim("forR_plotting/overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF-only")
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF.bound <- read.delim("forR_plotting/overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF-bound")
#Add pseudocount of 0.1 to number of overlaps
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_dELS$Num_overlaps <- overlapCounts_TEsubfamilies_fullClassCCREs_wClass_dELS$Num_overlaps + 0.1
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_pELS$Num_overlaps <- overlapCounts_TEsubfamilies_fullClassCCREs_wClass_pELS$Num_overlaps + 0.1
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_PLS$Num_overlaps <- overlapCounts_TEsubfamilies_fullClassCCREs_wClass_PLS$Num_overlaps + 0.1
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase.H3K4me3$Num_overlaps <- overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase.H3K4me3$Num_overlaps + 0.1
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase.only$Num_overlaps <- overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase.only$Num_overlaps + 0.1
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF.only$Num_overlaps <- overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF.only$Num_overlaps + 0.1
overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF.bound$Num_overlaps <- overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF.bound$Num_overlaps + 0.1
#Function for plotting TE classes separately and then combining
plot_class_count_separate <- function(input_df) {
  scale_max <- max(input_df$Num_overlaps, na.rm = TRUE)
  scale_min <- min(input_df$Num_overlaps, na.rm = TRUE)
  p1 <- ggplot(input_df[input_df$Class=='DNA',], aes(Subfamily, Celltype, fill=Num_overlaps)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max), trans='log10') +
    labs(x="", y="", fill="number of\nelements", title="DNA") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p2 <- ggplot(input_df[input_df$Class=='LINE',], aes(Subfamily, Celltype, fill=Num_overlaps)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max), trans='log10') +
    labs(x="", y="", fill="number of\nelements", title="LINE") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p3 <- ggplot(input_df[input_df$Class=='SINE',], aes(Subfamily, Celltype, fill=Num_overlaps)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max), trans='log10') +
    labs(x="", y="", fill="number of\nelements", title="SINE") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p4 <- ggplot(input_df[input_df$Class=='LTR',], aes(Subfamily, Celltype, fill=Num_overlaps)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max), trans='log10') +
    labs(x="", y="", fill="number of\nelements", title="LTR") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  p5 <- ggplot(input_df[input_df$Class=='SVA',], aes(Subfamily, Celltype, fill=Num_overlaps)) +
    geom_tile() +
    scale_fill_gradient2(low='blue', mid="white", high="red", na.value='grey', limits=c(scale_min, scale_max), trans='log10') +
    labs(x="", y="", fill="number of\nelements", title="SVA") +
    theme(legend.title.align = 0.5, axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  lay <- rbind(c(1, 1, 1, 1, NA, NA), c(2, 2, 3, 5, NA, NA), c(4, 4, 4, 4, 4, 4))
  grid.arrange(p1, p2, p3, p4, p5, ncol=1, layout_matrix=lay)
}
#Plot (save as 30x100?)
plot_class_count_separate(overlapCounts_TEsubfamilies_fullClassCCREs_wClass_dELS)
plot_class_count_separate(overlapCounts_TEsubfamilies_fullClassCCREs_wClass_pELS)
plot_class_count_separate(overlapCounts_TEsubfamilies_fullClassCCREs_wClass_PLS)
plot_class_count_separate(overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase.H3K4me3)
plot_class_count_separate(overlapCounts_TEsubfamilies_fullClassCCREs_wClass_DNase.only)
plot_class_count_separate(overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF.only)
plot_class_count_separate(overlapCounts_TEsubfamilies_fullClassCCREs_wClass_CTCF.bound)

