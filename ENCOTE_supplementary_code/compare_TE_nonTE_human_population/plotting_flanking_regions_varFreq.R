#R 3.5.1
library(ggplot2)
source('../reference_files/ggplotting_source', encoding = 'UTF-8')

#Load empirical p-values
empirical_pvals_flanks_nonZero_shuffled_mean_diff <- read.delim("empirical_pvals_flanks_nonZero_shuffled_mean_diff")
empirical_pvals_flanks_nonZero_shuffled_nonTE_mean <- read.delim("empirical_pvals_flanks_nonZero_shuffled_nonTE_mean")
empirical_pvals_flanks_nonZero_shuffled_nonTE_perc <- read.delim("empirical_pvals_flanks_nonZero_shuffled_nonTE_perc")
empirical_pvals_flanks_nonZero_shuffled_perc_diff <- read.delim("empirical_pvals_flanks_nonZero_shuffled_perc_diff")
empirical_pvals_flanks_nonZero_shuffled_TE_mean <- read.delim("empirical_pvals_flanks_nonZero_shuffled_TE_mean")
empirical_pvals_flanks_nonZero_shuffled_TE_perc <- read.delim("empirical_pvals_flanks_nonZero_shuffled_TE_perc")
empirical_pvals_flanks_nonZero_results <- read.delim("empirical_pvals_flanks_nonZero_results")
#Plot (save as 4x17)
#Mean difference between TE and non-TE cCREs
observed_mean_diff <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='mean_diff')
ggplot(empirical_pvals_flanks_nonZero_shuffled_mean_diff, aes(Mean_diff)) +
  geom_histogram(binwidth=0.001) +
  labs(x='Difference in means', y='Shuffle count') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(data=observed_mean_diff, aes(xintercept=Observed), color='red') +
  facet_wrap(~cCRE_type, nrow=1)
#Percent overlap difference between TE and non-TE cCREs
observed_perc_diff <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='perc_diff')
ggplot(empirical_pvals_flanks_nonZero_shuffled_perc_diff, aes(Perc_diff)) +
  geom_histogram(binwidth=0.001) +
  labs(x='Difference in percentage', y='Shuffle count') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(data=observed_perc_diff, aes(xintercept=Observed), color='red') +
  facet_wrap(~cCRE_type, nrow=1)
#Mean difference between TE-cCREs and flanks
observed_TE_mean_diff_left <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='TE_mean_left')
observed_TE_mean_diff_left$Shuffle_dir <- 'left'
observed_TE_mean_diff_right <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='TE_mean_right')
observed_TE_mean_diff_right$Shuffle_dir <- 'right'
observed_TE_mean_diff <- rbind(observed_TE_mean_diff_left, observed_TE_mean_diff_right)
ggplot(empirical_pvals_flanks_nonZero_shuffled_TE_mean, aes(Difference)) +
  geom_histogram(binwidth=0.001) +
  labs(x='Difference in means', y='Shuffle count') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(data=observed_TE_mean_diff, aes(xintercept=Observed), color='red') +
  facet_grid(cols=vars(cCRE_type), rows=vars(Shuffle_dir))
#Percent overlap difference between TE-cCREs and flanks
observed_TE_perc_diff_left <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='TE_perc_left')
observed_TE_perc_diff_left$Shuffle_dir <- 'left'
observed_TE_perc_diff_right <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='TE_perc_right')
observed_TE_perc_diff_right$Shuffle_dir <- 'right'
observed_TE_perc_diff <- rbind(observed_TE_perc_diff_left, observed_TE_perc_diff_right)
ggplot(empirical_pvals_flanks_nonZero_shuffled_TE_perc, aes(Difference)) +
  geom_histogram(binwidth=0.001) +
  labs(x='Difference in percentage with variant overlap', y='Shuffle count') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(data=observed_TE_perc_diff, aes(xintercept=Observed), color='red') +
  facet_grid(cols=vars(cCRE_type), rows=vars(Shuffle_dir))
#Mean difference between non-TE cCREs and flanks
observed_nonTE_mean_diff_left <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='nonTE_mean_left')
observed_nonTE_mean_diff_left$Shuffle_dir <- 'left'
observed_nonTE_mean_diff_right <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='nonTE_mean_right')
observed_nonTE_mean_diff_right$Shuffle_dir <- 'right'
observed_nonTE_mean_diff <- rbind(observed_nonTE_mean_diff_left, observed_nonTE_mean_diff_right)
ggplot(empirical_pvals_flanks_nonZero_shuffled_nonTE_mean, aes(Difference)) +
  geom_histogram(binwidth=0.001) +
  labs(x='Difference in means', y='Shuffle count') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(data=observed_nonTE_mean_diff, aes(xintercept=Observed), color='red') +
  facet_grid(cols=vars(cCRE_type), rows=vars(Shuffle_dir))
#Percent overlap difference between non-TE cCREs and flanks
observed_nonTE_perc_diff_left <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='nonTE_perc_left')
observed_nonTE_perc_diff_left$Shuffle_dir <- 'left'
observed_nonTE_perc_diff_right <- subset(empirical_pvals_flanks_nonZero_results, Comparison=='nonTE_perc_right')
observed_nonTE_perc_diff_right$Shuffle_dir <- 'right'
observed_nonTE_perc_diff <- rbind(observed_nonTE_perc_diff_left, observed_nonTE_perc_diff_right)
ggplot(empirical_pvals_flanks_nonZero_shuffled_nonTE_perc, aes(Difference)) +
  geom_histogram(binwidth=0.001) +
  labs(x='Difference in percentage with variant overlap', y='Shuffle count') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(data=observed_nonTE_perc_diff, aes(xintercept=Observed), color='red') +
  facet_grid(cols=vars(cCRE_type), rows=vars(Shuffle_dir))

#Percent overlapping variants
flank_varFreq_forOverlapPercent <- read.delim("forR_plotting/flank_varFreq_forOverlapPercent")
#Save as 4x8.5
ggplot(flank_varFreq_forOverlapPercent, aes(Annotation, Percent, fill=Annotation)) +
  stat_summary(fun.data = mean_sdl, geom='bar') +
  geom_jitter(width=0, height=0) +
  scale_fill_manual(values=c('flank' = 'gray', 'non-TE' = '#F8766D', 'TE' = '#00BFC4')) +
  coord_cartesian(ylim=c(60, 72)) +
  labs(x='', y='Percent overlapping variants') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='none') +
  facet_wrap(~cCRE_type, nrow=1)

#Variant frequency
flank_varFreq_forDistribution <- read.delim("forR_plotting/flank_varFreq_forDistribution")
#Save as
ggplot(flank_varFreq_forDistribution, aes(Annotation, varFreq, fill=Annotation)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", color='black') +
  scale_fill_manual(values=c('flank' = 'gray', 'non-TE' = '#F8766D', 'TE' = '#00BFC4')) +
  labs(x='', y='Variants per 100bp') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='none') +
  facet_wrap(~cCRE_type, nrow=1)
ggplot(flank_varFreq_forDistribution, aes(Annotation, varFreq, fill=Annotation)) +
  geom_boxplot(outlier.shape=NA) +
  stat_summary(fun=mean, geom="point", color='black') +
  scale_fill_manual(values=c('flank' = 'gray', 'non-TE' = '#F8766D', 'TE' = '#00BFC4')) +
  coord_cartesian(ylim=c(0,2)) +
  labs(x='', y='Variants per 100bp') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='none') +
  facet_wrap(~cCRE_type, nrow=1)
ggplot(flank_varFreq_forDistribution, aes(varFreq, color=Annotation)) +
  geom_freqpoly(binwidth=0.1) +
  scale_color_manual(values=c('flank' = 'gray', 'non-TE' = '#F8766D', 'TE' = '#00BFC4')) +
  scale_y_continuous(trans='log2') +
  labs(x='Variants per 100bp') +
  facet_wrap(~cCRE_type)

#Means
mean(subset(flank_varFreq_forDistribution, Annotation == 'TE' & cCRE_type == 'PLS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'non-TE' & cCRE_type == 'PLS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'flank' & cCRE_type == 'PLS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'TE' & cCRE_type == 'pELS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'non-TE' & cCRE_type == 'pELS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'flank' & cCRE_type == 'pELS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'TE' & cCRE_type == 'dELS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'non-TE' & cCRE_type == 'dELS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'flank' & cCRE_type == 'dELS')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'TE' & cCRE_type == 'DNase-H3K4me3')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'non-TE' & cCRE_type == 'DNase-H3K4me3')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'flank' & cCRE_type == 'DNase-H3K4me3')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'TE' & cCRE_type == 'CTCF-only')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'non-TE' & cCRE_type == 'CTCF-only')$varFreq)
mean(subset(flank_varFreq_forDistribution, Annotation == 'flank' & cCRE_type == 'CTCF-only')$varFreq)
