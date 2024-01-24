#R 3.5.1
library(ggplot2)
library(dplyr)
source('reference_files/ggplotting_source', encoding = 'UTF-8')

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#TE proportions
proportion_genomic_TEs <- read.delim("proportion_genomic_TEs")
proportion_genomic_TEs$Class <- factor(proportion_genomic_TEs$Class,
                                       levels=c("DNA", "LINE", "LTR", "SINE", "SVA", "non-TE"))
#Compute the position of labels
data <- proportion_genomic_TEs %>% 
  arrange(desc(Class)) %>%
  mutate(prop = Genomic_proportion / sum(proportion_genomic_TEs$Genomic_proportion) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop )
data$xpos <- 1
#Plot piechart
ggplot(data, aes(x="", y=prop, fill=Class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "SVA"="#E76BF3", "non-TE"="grey")) +
  geom_text(aes(x= xpos, y = ypos, label = Genomic_percentage), color = "black", size=6) +
  theme_void()
#Plot piechart without SVA
data2 <- subset(data, Class!="SVA")
ggplot(data2, aes(x="", y=prop, fill=Class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "non-TE"="grey")) +
  geom_text(aes(x= xpos, y = ypos, label = Genomic_percentage), color = "black", size=6) +
  theme_void()

#TE-cCRE overlap
proportion_cCREs_classes_justTEs <- read.delim("proportion_cCREs_classes_justTEs", comment.char="#")
proportion_cCREs_classes_justTEs$Class <- factor(proportion_cCREs_classes_justTEs$Class,
                                                 levels=c("non-TE", "TE",
                                                          "DNA", "LINE", "LTR", "SINE", "SVA"))
proportion_cCREs_classes_justTEs <- subset(proportion_cCREs_classes_justTEs, Class != "TE")
#Add genomic TE proportions as background/expectation
colnames(proportion_genomic_TEs)[2] <- "cCRE_num_perc"
proportion_genomic_TEs$cCRE_type <- "Genomic"
proportion_genomic_TEs$cCRE_num <- "NA"
proportion_genomic_TEs$cCRE_bp <- "NA"
proportion_genomic_TEs$cCRE_bp_perc <- proportion_genomic_TEs$cCRE_num_perc
proportion_cCREs_classes_justTEs$Genomic_percentage <- "NA"
proportion_cCREs_classes_justTEs <- rbind(proportion_cCREs_classes_justTEs, proportion_genomic_TEs)
proportion_cCREs_classes_justTEs$cCRE_type <- factor(proportion_cCREs_classes_justTEs$cCRE_type,
                                                 levels=c("Genomic", "All cCREs", "dELS", "pELS",
                                                          "PLS", "DNase-H3K4me3", "CTCF-only", "CTCF-bound"))
#Plot with SVA
ggplot(proportion_cCREs_classes_justTEs, aes(x=cCRE_type, y=cCRE_num_perc, fill=Class)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "SVA"="#E76BF3",
                                         "non-TE"="grey")) +
  #scale_fill_manual("TE class", values=c("grey", "orange", "white", "white", "white", "white",
  #                                       "#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),
  #                  drop=FALSE) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')
#Plot without SVA
proportion_cCREs_classes_justTEs_woSVA <- subset(proportion_cCREs_classes_justTEs, Class != "SVA")
ggplot(proportion_cCREs_classes_justTEs_woSVA, aes(x=cCRE_type, y=cCRE_num_perc, fill=Class)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", 
                                         "non-TE"="grey")) +
  #scale_fill_manual("TE class", values=c("grey", "orange", "white", "white", "white", "white",
  #                                       "#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),
  #                  drop=FALSE) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')

#Plots without CTCF-bound group
#Plot with SVA
proportion_cCREs_classes_justTEs_noBound <- subset(proportion_cCREs_classes_justTEs, cCRE_type != "CTCF-bound")
ggplot(proportion_cCREs_classes_justTEs_noBound, aes(x=cCRE_type, y=cCRE_num_perc, fill=Class)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "SVA"="#E76BF3",
                                         "non-TE"="grey")) +
  #scale_fill_manual("TE class", values=c("grey", "orange", "white", "white", "white", "white",
  #                                       "#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),
  #                  drop=FALSE) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')
#Plot without SVA
proportion_cCREs_classes_justTEs_woSVA_noBound <- subset(proportion_cCREs_classes_justTEs_noBound, Class != "SVA")
ggplot(proportion_cCREs_classes_justTEs_woSVA_noBound, aes(x=cCRE_type, y=cCRE_num_perc, fill=Class)) +
  geom_bar(position="fill", stat="identity", width=0.5) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", 
                                         "non-TE"="grey")) +
  #scale_fill_manual("TE class", values=c("grey", "orange", "white", "white", "white", "white",
  #                                       "#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),
  #                  drop=FALSE) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')
#Plot without non-TE as well (saved as 6x12 for retreat)
proportion_cCREs_classes_justTEs_woSVA_noBound_woNonTE <- subset(proportion_cCREs_classes_justTEs_woSVA_noBound, Class != "non-TE")
ggplot(proportion_cCREs_classes_justTEs_woSVA_noBound_woNonTE, aes(x=cCRE_type, y=cCRE_num_perc, fill=Class)) +
  geom_bar(stat="identity", width=0.5) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6")) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')

#Full classification cell types TE-cCRE overlap (by class, split by cCRE type)
#First do all cCREs
proportions_TEclass_fullClassCCREs_All_cCREs <- read.delim("forR_plotting/proportions_TEclass_fullClassCCREs_All_cCREs")
proportions_TEclass_fullClassCCREs_All_cCREs <- subset(proportions_TEclass_fullClassCCREs_All_cCREs,
                                                       Class != 'Other' & Class != 'All_REs' & Class != 'TE')
proportions_TEclass_fullClassCCREs_All_cCREs$Class <- factor(proportions_TEclass_fullClassCCREs_All_cCREs$Class,
                                       levels=c("non-TE", "DNA", "LINE", "LTR", "SINE", "SVA"))
ggplot(proportions_TEclass_fullClassCCREs_All_cCREs, aes(Cell_type, cCRE_num_perc, fill=Class)) +
  geom_bar(position="fill", stat="identity", width=1) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "SVA"="#E76BF3",
                                         "non-TE"="grey")) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')
#Plot without SVA and non-TE
proportions_TEclass_fullClassCCREs_All_cCREs_woSVA_woNonTE <- subset(proportions_TEclass_fullClassCCREs_All_cCREs, Class != "SVA" & Class != "non-TE")
ggplot(proportions_TEclass_fullClassCCREs_All_cCREs_woSVA_woNonTE, aes(Cell_type, cCRE_num_perc, fill=Class)) +
  geom_bar(stat="identity", width=0.5) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6")) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')
#Simplify cell/tissue name for Full classification cell types TE-cCRE overlap (by class, split by cCRE type)
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange <- read.delim("forR_plotting/proportions_TEclass_fullClassCCREs_All_cCREs_nameChange")
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange <- subset(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange,
                                                       Class != 'Other' & Class != 'All_REs' & Class != 'TE')
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange$Class <- factor(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange$Class,
                                                             levels=c("non-TE", "DNA", "LINE", "LTR", "SINE", "SVA"))
cell_types_sorted <- c(
  #Primary cells
  'astrocyte', 'B cell', 'CD14-positive monocyte', 'fibroblast of dermis', 'keratinocyte', 
  #In vitro differentiated cells
  'bipolar neuron', 'hepatocyte (embryonic)','myotube', 'neural progenitor cell', 
  #Cell lines (embryonic/normal)
  'AG04450', 'GM12878', 'GM23338', 'H1 hESC', 'IMR-90', 
  #Cell lines (cancer)
  'A673', 'HCT116', 'HeLa-S3', 'HepG2', 'K562', 'MCF-7', 'MM.1S', 'OCI-LY7', 'Panc1', 'PC-3', 'PC-9'
                       )
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange$Cell_type <- factor(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange$Cell_type,
                                                                            levels=cell_types_sorted)
#Remove SVA and non-TE
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_woSVA_woNonTE <- subset(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange, Class != "SVA" & Class != "non-TE")
ggplot(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_woSVA_woNonTE, aes(Cell_type, cCRE_num_perc, fill=Class)) +
  geom_bar(stat="identity", width=0.5) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6")) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')
#Plot TE proportions only
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo <- read.delim("forR_plotting/proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo")
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo <- subset(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo,
                                                                  Class != 'Other' & Class != 'All_REs' & Class != 'TE')
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo$Class <- factor(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo$Class,
                                                                        levels=c("non-TE", "DNA", "LINE", "LTR", "SINE", "SVA"))
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo$Cell_type <- factor(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo$Cell_type,
                                                                            levels=cell_types_sorted)
proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo <- subset(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo, Class != "SVA" & Class != "non-TE")
ggplot(proportions_TEclass_fullClassCCREs_All_cCREs_nameChange_moreInfo, aes(Cell_type, TE_cCRE_num_perc, fill=Class)) +
  geom_col(width=0.5) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6")) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y='Annotated proportion', x='')

