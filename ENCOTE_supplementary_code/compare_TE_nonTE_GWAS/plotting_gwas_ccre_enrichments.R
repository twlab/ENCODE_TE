#R 3.5.1
library(ggplot2)
source('../reference_files/ggplotting_source', encoding = 'UTF-8')

#Load data
enrichments_gwas_cCRE_overlap <- read.delim("enrichments_gwas_cCRE_overlap")
#Shuffled count distributions
ggplot(enrichments_gwas_cCRE_overlap, aes(Shuffled_count)) +
  geom_density() +
  facet_grid(rows=vars(Parent_term), cols=vars(cCRE_group), scales='free')
#Enrichment distributions
ggplot(enrichments_gwas_cCRE_overlap, aes(Parent_term, Enrichment, fill=cCRE_group)) +
  geom_boxplot() +
  geom_hline(yintercept=1, color='grey', linetype='dashed') +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'top') +
  labs(x='', y='Enrichment')
#Select parent terms
enrichments_gwas_cCRE_overlap_selectTerms <- subset(enrichments_gwas_cCRE_overlap, Parent_term != 'Biological process' &
                                                      Parent_term != 'NR' & Parent_term != 'Other disease' &
                                                      Parent_term != 'Other measurement' & Parent_term != 'Other trait')
ggplot(enrichments_gwas_cCRE_overlap_selectTerms, aes(Parent_term, Enrichment, fill=cCRE_group)) +
  geom_boxplot() +
  geom_hline(yintercept=1, color='grey', linetype='dashed') +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'top') +
  labs(x='', y='Enrichment')
#Select groupings (All cCRE, non-TE, and TE)
enrichments_gwas_cCRE_overlap_selectTerms_selectGroups <- subset(enrichments_gwas_cCRE_overlap_selectTerms, 
                                                                 cCRE_group != 'human_nonTE' & cCRE_group != 'human_TE')
ggplot(enrichments_gwas_cCRE_overlap_selectTerms_selectGroups, aes(Parent_term, Enrichment, fill=cCRE_group)) +
  geom_boxplot() +
  geom_hline(yintercept=1, color='grey', linetype='dashed') +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'top') +
  labs(x='', y='Enrichment')

#Negative control (shuffled)
negative_shuffled_enrichments <- read.delim("negative_shuffled_enrichments")
ggplot(negative_shuffled_enrichments, aes(Parent_term, Enrichment, fill=cCRE_group)) +
  geom_boxplot() +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'top') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'top') +
  geom_hline(yintercept=1, color='grey', linetype='dashed') +
  labs(x='', y='Enrichment')
