#R 3.5.1
library(ggplot2)
source('reference_files/ggplotting_source', encoding = 'UTF-8')

#Load data
final_info_human <- read.delim2("final_info_human")
orthologousTE_info <- subset(final_info_human, ! is.na(Human_TE) & Synteny_annotation == 'TEortholog')
orthologousTE_info$Annotation <- "Orthologous TE"
TE_info <- subset(final_info_human, ! is.na(Human_TE)  & Synteny_annotation != 'TEortholog')
TE_info$Annotation <- "TE Human-lineage"
syntenic_info <- subset(final_info_human, is.na(Human_TE) & Synteny_annotation == 'Syntenic')
humanOnly_info <- subset(final_info_human, Synteny_annotation == 'Human_specific')

#Compare kimura divergence (age) between orthologous TEs and non-orthologous TEs
compare_orthologous_vs_non <- rbind(orthologousTE_info, TE_info)
compare_orthologous_vs_non$Human_TE_age <- as.numeric(as.character(compare_orthologous_vs_non$Human_TE_age))
ggplot(compare_orthologous_vs_non, aes(Human_TE_age, fill=Annotation)) +
  geom_density(alpha=0.5) +
  labs(x='Human TE kimura divergence', y='Density')

#Add clade info (only has human cCREs that are associated with human TE)
#Load data
final_info_human_wClade <- read.delim("final_info_human_wClade")
orthologousTE_info_wClade <- subset(final_info_human_wClade, Synteny_annotation == 'TEortholog')
orthologousTE_info_wClade$Annotation <- "Orthologous TE"
TE_info_wClade <- subset(final_info_human_wClade, Synteny_annotation != 'TEortholog')
TE_info_wClade$Annotation <- "Human-lineage TE"
#Compare kimura divergence (age) between orthologous TEs and non-orthologous TEs
compare_orthologous_vs_non <- rbind(orthologousTE_info_wClade, TE_info_wClade)
compare_orthologous_vs_non$Human_TE_age <- as.numeric(as.character(compare_orthologous_vs_non$Human_TE_age))
compare_orthologous_vs_non$Clade <- factor(compare_orthologous_vs_non$Clade, 
                                           levels=c('Metazoa', 'Vertebrata', 'Euteleostomi', 'Tetrapoda', 'Amniota', 
                                                    'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires',
                                                    'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 
                                                    'Hominoidea', 'Hominidae', 'Homininae', 'Homo sapiens'))
compare_orthologous_vs_non$Species_classification <- factor(compare_orthologous_vs_non$Species_classification,
                                                            levels=c('Shared_subfamily', 'Lineage_subfamily'))
ggplot(compare_orthologous_vs_non, aes(Human_TE_age, fill=Annotation)) +
  geom_density(alpha=0.5) +
  labs(x='Human TE kimura divergence %', y='Density')
compare_orthologous_vs_non_noNAclassifications <- subset(compare_orthologous_vs_non, ! is.na(Species_classification))
ggplot(compare_orthologous_vs_non_noNAclassifications, aes(x=Annotation, fill=Species_classification)) +
  geom_bar(position='fill') +
  labs(x='', y='Percentage')
ggplot(compare_orthologous_vs_non_noNAclassifications, aes(x=Annotation, fill=Clade)) +
  geom_bar(position='fill') +
  labs(x='', y='Percentage')

#Add phyloP and phastCons scores
final_info_human_30way <- read.delim("phylop_phastcons/forR_final_info_human_30way")
final_info_human_100way <- read.delim("phylop_phastcons/forR_final_info_human_100way")
#Plot
ggplot(final_info_human_100way, aes(Human_cCRE_type, phyloP_score)) +
  geom_violin() +
  labs(x='', y='phyloP100 way score') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(Synteny_annotation ~ cCRE_comparison)
ggplot(final_info_human_100way, aes(Human_cCRE_type, phastCons_score)) +
  geom_violin() +
  labs(x='', y='phastCons100 way score') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(Synteny_annotation ~ cCRE_comparison)
ggplot(final_info_human_30way, aes(Human_cCRE_type, phyloP_score)) +
  geom_violin() +
  labs(x='', y='phyloP30 way score') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(Synteny_annotation ~ cCRE_comparison)
ggplot(final_info_human_30way, aes(Human_cCRE_type, phastCons_score)) +
  geom_violin() +
  labs(x='', y='phastCons30 way score') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(Synteny_annotation ~ cCRE_comparison)

##### Conservation by phyloP and phastCons, comparing orthologous TE vs. non-orthologous #####
#phyloP and phastCons, TEs only (with clade info)
final_info_human_wClade_30way <- read.delim("phylop_phastcons/forR_final_info_human_wClade_30way")
final_info_human_wClade_100way <- read.delim("phylop_phastcons/forR_final_info_human_wClade_100way")
#Split into orthologous TE and human-lineage TE
orthologousTE_info_wClade_30way <- subset(final_info_human_wClade_30way, Synteny_annotation == 'TEortholog')
orthologousTE_info_wClade_30way$Annotation <- "Orthologous TE"
TE_info_wClade_30way <- subset(final_info_human_wClade_30way, Synteny_annotation != 'TEortholog')
TE_info_wClade_30way$Annotation <- "Human-lineage TE"
orthologousTE_info_wClade_100way <- subset(final_info_human_wClade_100way, Synteny_annotation == 'TEortholog')
orthologousTE_info_wClade_100way$Annotation <- "Orthologous TE"
TE_info_wClade_100way <- subset(final_info_human_wClade_100way, Synteny_annotation != 'TEortholog')
TE_info_wClade_100way$Annotation <- "Human-lineage TE"
compare_orthologous_vs_non_30way <- rbind(orthologousTE_info_wClade_30way, TE_info_wClade_30way)
compare_orthologous_vs_non_30way$Human_TE_age <- as.numeric(as.character(compare_orthologous_vs_non_30way$Human_TE_age))
compare_orthologous_vs_non_30way$Clade <- factor(compare_orthologous_vs_non_30way$Clade, 
                                           levels=c('Metazoa', 'Vertebrata', 'Euteleostomi', 'Tetrapoda', 'Amniota', 
                                                    'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires',
                                                    'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 
                                                    'Hominoidea', 'Hominidae', 'Homininae', 'Homo sapiens'))
compare_orthologous_vs_non_30way$Species_classification <- factor(compare_orthologous_vs_non_30way$Species_classification,
                                                            levels=c('Shared_subfamily', 'Lineage_subfamily'))
compare_orthologous_vs_non_100way <- rbind(orthologousTE_info_wClade_100way, TE_info_wClade_100way)
compare_orthologous_vs_non_100way$Human_TE_age <- as.numeric(as.character(compare_orthologous_vs_non_100way$Human_TE_age))
compare_orthologous_vs_non_100way$Clade <- factor(compare_orthologous_vs_non_100way$Clade, 
                                                 levels=c('Metazoa', 'Vertebrata', 'Euteleostomi', 'Tetrapoda', 'Amniota', 
                                                          'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires',
                                                          'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 
                                                          'Hominoidea', 'Hominidae', 'Homininae', 'Homo sapiens'))
compare_orthologous_vs_non_100way$Species_classification <- factor(compare_orthologous_vs_non_100way$Species_classification,
                                                                  levels=c('Shared_subfamily', 'Lineage_subfamily'))
compare_orthologous_vs_non_noNAclassifications_30way <- subset(compare_orthologous_vs_non_30way, ! is.na(Species_classification))
compare_orthologous_vs_non_noNAclassifications_100way <- subset(compare_orthologous_vs_non_100way, ! is.na(Species_classification))
#Plot
#Compare shared vs. lineage specific subfamilies
ggplot(compare_orthologous_vs_non_noNAclassifications_100way, aes(Synteny_annotation, phyloP_score, fill=Species_classification)) +
  geom_violin() +
  labs(x='', y='phyloP100 way score') +
  facet_grid(rows=vars(Human_cCRE_type))
ggplot(compare_orthologous_vs_non_noNAclassifications_100way, aes(Synteny_annotation, phastCons_score, fill=Species_classification)) +
  geom_violin() +
  labs(x='', y='phastCons100 way score') +
  facet_grid(rows=vars(Human_cCRE_type))
ggplot(compare_orthologous_vs_non_noNAclassifications_30way, aes(Synteny_annotation, phyloP_score, fill=Species_classification)) +
  geom_violin() +
  labs(x='', y='phyloP30 way score') +
  facet_grid(rows=vars(Human_cCRE_type))
ggplot(compare_orthologous_vs_non_noNAclassifications_30way, aes(Synteny_annotation, phastCons_score, fill=Species_classification)) +
  geom_violin() +
  labs(x='', y='phastCons way score') +
  facet_grid(rows=vars(Human_cCRE_type))
#Compare 

##### Compare orthologous TE vs. non-TE shared cCRE conservation #####
#Make sure final_info_human_30way and final_info_human_100way are loaded
compare_TE_vs_non_30way <- subset(final_info_human_30way, cCRE_comparison == 'Shared_cCRE')
compare_TE_vs_non_100way <- subset(final_info_human_100way, cCRE_comparison == 'Shared_cCRE')
ggplot(compare_TE_vs_non_30way, aes(Synteny_annotation, phastCons_score)) +
  geom_violin() +
  labs(x='', y='phastCons30 way score') +
  facet_grid(rows=vars(Human_cCRE_type))
ggplot(compare_TE_vs_non_100way, aes(Synteny_annotation, phastCons_score)) +
  geom_violin() +
  labs(x='', y='phastCons100 way score') +
  facet_grid(rows=vars(Human_cCRE_type))

