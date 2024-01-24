#R 3.5.1
library(ggplot2)
library(gridExtra)
source('reference_files/ggplotting_source', encoding = 'UTF-8')

#Load data
distribution_TEclass_subfamily_celltype_enrichments <- read.delim("forR_plotting/distribution_TEclass_subfamily_celltype_enrichments")

#ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_colors <- gg_color_hue(6)

#Plot each class separately, getting rid of 0 cell type overlap
ggplot(distribution_TEclass_subfamily_celltype_enrichments, aes(NumEnrichedCelltypes, NumSubfamilies, fill=Class)) +
  geom_col() +
  xlim(c(1,NA)) +
  labs(x='Number of cell types (enriched)', y='Number of subfamilies', title='Distribution of enriched subfamilies by class') +
  facet_wrap(Class~cCRE, scales='free_y', nrow=5)
#Plot each class separately
p1 <- ggplot(distribution_TEclass_subfamily_celltype_enrichments[distribution_TEclass_subfamily_celltype_enrichments$Class=='DNA',], aes(NumEnrichedCelltypes, NumSubfamilies)) +
  geom_col(fill=gg_colors[1]) +
  labs(x='', y='', title='DNA') +
  #scale_x_continuous(limits=c(1,NA)) +
  facet_wrap(~cCRE, nrow=1)
p2 <- ggplot(distribution_TEclass_subfamily_celltype_enrichments[distribution_TEclass_subfamily_celltype_enrichments$Class=='LINE',], aes(NumEnrichedCelltypes, NumSubfamilies)) +
  geom_col(fill=gg_colors[2]) +
  labs(x='', y='', title='LINE') +
  #scale_x_continuous(limits=c(1,NA)) +
  facet_wrap(~cCRE, nrow=1)
p3 <- ggplot(distribution_TEclass_subfamily_celltype_enrichments[distribution_TEclass_subfamily_celltype_enrichments$Class=='LTR',], aes(NumEnrichedCelltypes, NumSubfamilies)) +
  geom_col(fill=gg_colors[3]) +
  labs(x='', y='Number of subfamilies', title='LTR') +
  #scale_x_continuous(limits=c(1,NA)) +
  facet_wrap(~cCRE, nrow=1)
p4 <- ggplot(distribution_TEclass_subfamily_celltype_enrichments[distribution_TEclass_subfamily_celltype_enrichments$Class=='SINE',], aes(NumEnrichedCelltypes, NumSubfamilies)) +
  geom_col(fill=gg_colors[4]) +
  labs(x='', y='', title='SINE') +
  #scale_x_continuous(limits=c(1,NA)) +
  facet_wrap(~cCRE, nrow=1)
p5 <- ggplot(distribution_TEclass_subfamily_celltype_enrichments[distribution_TEclass_subfamily_celltype_enrichments$Class=='SVA',], aes(NumEnrichedCelltypes, NumSubfamilies)) +
  geom_col(fill=gg_colors[5]) +
  labs(x='Number of cell types (enriched)', y='', title='SVA') +
  #scale_x_continuous(limits=c(1,NA)) +
  facet_wrap(~cCRE, nrow=1)
grid.arrange(p1, p2, p3, p4, p5, ncol=1)


