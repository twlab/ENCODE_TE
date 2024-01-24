#R 3.5.1
library(ggplot2)
library(ggpubr)
library(reshape2)
source('~/ggplotting_source', encoding = 'UTF-8')

#Randomly selected background, enriched motif ancestral origin
#Load data
regionOnly_originAnc_confidentMotifs_consensus_perSubfamily_sepInts_compareRandom <- read.delim("final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perSubfamily_sepInts_compareRandom")
#Plot
observed_data_filtered <- subset(regionOnly_originAnc_confidentMotifs_consensus_perSubfamily_sepInts_compareRandom,
                                 Class != "SVA" & Class != "Other")
observed_data_filtered$Class <- factor(observed_data_filtered$Class, levels=c("DNA", "LINE", "LTR", "ERV-int", "SINE"))
ggplot(observed_data_filtered, aes(Class, Perc_anc_origin, fill=Class)) +
  geom_violin(bw=10) +
  geom_boxplot(width=0.2) +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  labs(x='', y='Percent ancestral origin', title='Consensus-based origin of confident motifs in cCREs')
ggplot(observed_data_filtered, aes(Kimura_div, Perc_anc_origin, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Percent ancestral origin', title='Consensus-based origin of confident motifs in cCREs') +
  stat_cor(label.y = 100, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), vjust=-0.5) + 
  stat_regline_equation(label.y = 100) +
  facet_wrap(~Class, nrow=1)
#Save as 4x16
ggplot(observed_data_filtered, aes(Kimura_div, Perc_anc_origin, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Percent ancestral origin') +
  stat_cor(label.y = 95, label.x=20, aes(label = ..rr.label..), vjust=-0.5) + 
  stat_cor(label.y = 95, label.x=20, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class, nrow=1)

summary(lm(Perc_anc_origin ~ Kimura_div, data=observed_data_filtered, 
           subset=(Class=='LINE')))
summary(lm(Perc_anc_origin ~ Kimura_div, data=observed_data_filtered, 
           subset=(Class=='LTR')))
summary(lm(Perc_anc_origin ~ Kimura_div, data=observed_data_filtered, 
           subset=(Class=='ERV-int')))
summary(lm(Perc_anc_origin ~ Kimura_div, data=observed_data_filtered, 
           subset=(Class=='DNA')))
summary(lm(Perc_anc_origin ~ Kimura_div, data=observed_data_filtered, 
           subset=(Class=='SINE')))

library(dunn.test)
dunn.test(observed_data_filtered$Perc_anc_origin, observed_data_filtered$Class, method='bh')

#Number of enriched archetypes per subfamily
ggplot(observed_data_filtered, aes(Num_archetypes, fill=Class)) +
  geom_histogram(binwidth = 1)
ggplot(observed_data_filtered, aes(Num_archetypes, fill=Class)) +
  geom_histogram(binwidth=1) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,150)) +
  scale_x_continuous(breaks=seq(1,23))
#Number of enriched archetypes relation to age
ggplot(observed_data_filtered, aes(Kimura_div, Num_archetypes, color=Class)) +
  geom_point() +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Number of motif archetypes') +
  stat_cor(label.y = 17, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), vjust=-0.5) + 
  stat_regline_equation(label.y = 17, vjust=1) +
  facet_wrap(~Class, nrow=1)
#Number of enriched archetypes relation to length
length_data <- read.delim("final_ancOrigin_calculation_coverageControl/regionOnly_length_human_TEsubfamily_confidentMotifs_consensus_originAnc")
length_data_filtered <- subset(length_data, Class != "SVA" & Class != "Other")
length_data_filtered$Class <- factor(length_data_filtered$Class, levels=c("DNA", "LINE", "LTR", "ERV-int", "SINE"))
ggplot(length_data_filtered, aes(Length, Num_archetypes, color=Class)) +
  geom_point() +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  stat_smooth(method='lm') +
  labs(x='Length (bp)', y='Number of motif archetypes') +
  facet_wrap(~Class)

#####Randomly selected motifs (null expectation)#####
#Mean
ggplot(observed_data_filtered, aes(Kimura_div, Random_perc_anc, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Random percent ancestral origin') +
  stat_cor(label.y = 95, label.x=20, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), vjust=-0.5) + 
  stat_regline_equation(label.y = 95, label.x=20) +
  facet_wrap(~Class)
#Median
ggplot(observed_data_filtered, aes(Kimura_div, Median_random_perc_anc, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Random percent ancestral origin') +
  stat_cor(label.y = 95, label.x=20, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), vjust=-0.5) + 
  stat_regline_equation(label.y = 95, label.x=20) +
  facet_wrap(~Class)
#Compare observed cCRE enriched motifs to randomly selected motifs
temp_data <- melt(observed_data_filtered, id.vars=c('Subfamily', 'Class', 'Kimura_div', 'Num_archetypes'), variable.name='Perc_type')
temp_data$Perc_type <- sub('Perc_anc_origin', 'Observed', temp_data$Perc_type)
temp_data$Perc_type <- sub('Random_perc_anc', 'Random mean', temp_data$Perc_type)
temp_data$Perc_type <- sub('Median_random_perc_anc', 'Random median', temp_data$Perc_type)
#Mean and median violin+boxplot (6x15)
ggplot(temp_data, aes(Perc_type, value, fill=Class)) +
  geom_violin(bw=10) +
  geom_boxplot(width=0.2) +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  labs(x='', y='Percent ancestral origin') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~Class, nrow=1)
wilcox.test(subset(observed_data_filtered, Class=="LINE")$Perc_anc_origin, subset(observed_data_filtered, Class="LINE")$Random_perc_anc, conf.int=TRUE)
wilcox.test(subset(observed_data_filtered, Class=="SINE")$Perc_anc_origin, subset(observed_data_filtered, Class="SINE")$Random_perc_anc)
wilcox.test(subset(observed_data_filtered, Class=="DNA")$Perc_anc_origin, subset(observed_data_filtered, Class="DNA")$Random_perc_anc)
wilcox.test(subset(observed_data_filtered, Class=="LTR")$Perc_anc_origin, subset(observed_data_filtered, Class="LTR")$Random_perc_anc, conf.int=TRUE)
wilcox.test(subset(observed_data_filtered, Class=="ERV-int")$Perc_anc_origin, subset(observed_data_filtered, Class="ERV-int")$Random_perc_anc)
wilcox.test(subset(observed_data_filtered, Class=="LINE")$Perc_anc_origin, 
            subset(observed_data_filtered, Class="LINE")$Median_random_perc_anc, conf.int=TRUE)
wilcox.test(subset(observed_data_filtered, Class=="SINE")$Perc_anc_origin, 
            subset(observed_data_filtered, Class="SINE")$Median_random_perc_anc, conf.int=TRUE)
wilcox.test(subset(observed_data_filtered, Class=="DNA")$Perc_anc_origin, 
            subset(observed_data_filtered, Class="DNA")$Median_random_perc_anc, conf.int=TRUE)
wilcox.test(subset(observed_data_filtered, Class=="LTR")$Perc_anc_origin, 
            subset(observed_data_filtered, Class="LTR")$Median_random_perc_anc, conf.int=TRUE)
wilcox.test(subset(observed_data_filtered, Class=="ERV-int")$Perc_anc_origin, 
            subset(observed_data_filtered, Class="ERV-int")$Median_random_perc_anc, conf.int=TRUE)
p.adjust(c(5.227e-08, 0.8627, 4.683e-09, 2.2e-16, 0.005341), method='BH')
#Keep only median because it's more comparable to observation (a single value rather than mean of multiple values) (save as 6x12)
median_only_data <- subset(temp_data, Perc_type != "Random mean")
median_only_data$Perc_type <- sub('Random median', 'Random', median_only_data$Perc_type)
ggplot(median_only_data, aes(Perc_type, value, fill=Class)) +
  geom_violin(bw=10) +
  geom_boxplot(width=0.2) +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  labs(x='', y='Percent ancestral origin') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~Class, nrow=1)

##### Percentage of motifs in consensus sequences #####
motif_in_consensus <- read.delim("final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perMotif_nonConsensus")
motif_in_consensus$Class <- factor(motif_in_consensus$Class, levels=c("DNA", "LINE", "LTR", "ERV-int", "SINE"))
motif_in_consensus$Num_consensus_archetypes <- motif_in_consensus$Num_archetypes - motif_in_consensus$Num_nonConsensus_archetypes
motif_in_consensus$Perc_consensus <- motif_in_consensus$Num_consensus_archetypes/motif_in_consensus$Num_archetypes*100
motif_in_consensus <- subset(motif_in_consensus, Class != "SVA" & Class != "Other")
ggplot(motif_in_consensus, aes(Class, Perc_consensus, fill=Class)) +
  geom_violin(bw=10) +
  stat_summary(fun='median', geom='point', color='black') +
  #geom_boxplot(width=0.2) +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  labs(x='', y='Percent motifs found in consensus')
library(dunn.test)
dunn.test(motif_in_consensus$Perc_consensus, motif_in_consensus$Class, method='bh')


#Number of enriched archetypes per subfamily
ggplot(observed_data_filtered, aes(Num_archetypes, fill=Class)) +
  geom_histogram(binwidth=1) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,150)) +
  scale_x_continuous(breaks=seq(1,23))

