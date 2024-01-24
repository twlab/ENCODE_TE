#R 3.5.1
library(ggplot2)
library(ggpubr)
source('reference_files/ggplotting_source', encoding = 'UTF-8')

#####coverage and length control motifs (all) per Motif#####
#Load data
percent_motifs_lengthControl <- read.delim("final_ancOrigin_calculation_coverageControl/percent_motifs_lengthControl")
percent_motifs_coverageControl <- read.delim("final_ancOrigin_calculation_coverageControl/percent_motifs_coverageControl")
percent_motifs_lengthControl <- subset(percent_motifs_lengthControl,
                         Class != "SVA" & Class != "Other")
percent_motifs_lengthControl$Class <- factor(percent_motifs_lengthControl$Class, 
                                             levels=c("DNA", "LINE", "LTR", "ERV-int", "SINE"))
percent_motifs_coverageControl <- subset(percent_motifs_coverageControl,
                         Class != "SVA" & Class != "Other")
percent_motifs_coverageControl$Class <- factor(percent_motifs_coverageControl$Class, 
                                               levels=c("DNA", "LINE", "LTR", "ERV-int", "SINE"))
#Plot length controlled background percentages
ggplot(percent_motifs_lengthControl, aes(Kimura_div, Percent_motifs, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Percent elements with motif') +
  stat_cor(label.y = 95, label.x=20, aes(label = ..rr.label..), vjust=-0.5) + 
  stat_cor(label.y = 95, label.x=20, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class)
ggplot(percent_motifs_lengthControl, aes(Class, Percent_motifs, fill=Class)) +
  geom_violin(bw=10) +
  geom_boxplot(width=0.3) +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  labs(x='', y='Percent elements')
library(dunn.test)
dunn.test(percent_motifs_lengthControl$Percent_motifs, percent_motifs_lengthControl$Class, method='bh')
ggplot(percent_motifs_lengthControl, aes(Length, Percent_motifs, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  scale_x_log10() +
  stat_smooth(method='lm') +
  labs(x='Consensus length', y='Percent elements with motif') +
  #stat_cor(label.y = 95, label.x=1000, aes(label = ..rr.label..), vjust=-0.5) + 
  #stat_cor(label.y = 95, label.x=1000, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class)

#Plot coverage controlled background percentages
#Save as 4x16
ggplot(percent_motifs_coverageControl, aes(Kimura_div, Percent_motifs, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Percent elements with motif') +
  stat_cor(label.y = 95, label.x=20, aes(label = ..rr.label..), vjust=-0.5) + 
  stat_cor(label.y = 95, label.x=20, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class, nrow=1)
ggplot(percent_motifs_coverageControl, aes(Class, Percent_motifs, fill=Class)) +
  geom_violin(bw=10) +
  geom_boxplot(width=0.3) +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  labs(x='', y='Percent elements')
library(dunn.test)
dunn.test(percent_motifs_coverageControl$Percent_motifs, percent_motifs_coverageControl$Class, method='bh')
ggplot(percent_motifs_coverageControl, aes(Length, Percent_motifs, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  scale_x_log10() +
  stat_smooth(method='lm') +
  labs(x='Consensus length', y='Percent elements with motif') +
  #stat_cor(label.y = 95, label.x=1000, aes(label = ..rr.label..), vjust=-0.5) + 
  #stat_cor(label.y = 95, label.x=1000, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class)

#####coverage and length control motifs (all) per Subfamily#####
#Load data
percent_motifs_lengthControl_perSubfamily <- read.delim("final_ancOrigin_calculation_coverageControl/percent_motifs_lengthControl_perSubfamily")
percent_motifs_coverageControl_perSubfamily <- read.delim("final_ancOrigin_calculation_coverageControl/percent_motifs_coverageControl_perSubfamily")
percent_motifs_lengthControl_perSubfamily <- subset(percent_motifs_lengthControl_perSubfamily,
                                       Class != "SVA" & Class != "Other")
percent_motifs_lengthControl_perSubfamily$Class <- factor(percent_motifs_lengthControl_perSubfamily$Class, 
                                             levels=c("DNA", "LINE", "LTR", "ERV-int", "SINE"))
percent_motifs_coverageControl_perSubfamily <- subset(percent_motifs_coverageControl_perSubfamily,
                                         Class != "SVA" & Class != "Other")
percent_motifs_coverageControl_perSubfamily$Class <- factor(percent_motifs_coverageControl_perSubfamily$Class, 
                                               levels=c("DNA", "LINE", "LTR", "ERV-int", "SINE"))
#Plot length controlled background percentages
ggplot(percent_motifs_lengthControl_perSubfamily, aes(Kimura_div, Percent_motifs, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Percent elements with motif') +
  stat_cor(label.y = 95, label.x=20, aes(label = ..rr.label..), vjust=-0.5) + 
  stat_cor(label.y = 95, label.x=20, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class)
ggplot(percent_motifs_lengthControl_perSubfamily, aes(Class, Percent_motifs, fill=Class)) +
  geom_violin(bw=10) +
  geom_boxplot(width=0.3) +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  labs(x='', y='Percent elements')
library(dunn.test)
dunn.test(percent_motifs_lengthControl_perSubfamily$Percent_motifs, percent_motifs_lengthControl_perSubfamily$Class, method='bh')
ggplot(percent_motifs_lengthControl_perSubfamily, aes(Length, Percent_motifs, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  scale_x_log10() +
  stat_smooth(method='lm') +
  labs(x='Consensus length', y='Percent elements with motif') +
  #stat_cor(label.y = 95, label.x=1000, aes(label = ..rr.label..), vjust=-0.5) + 
  #stat_cor(label.y = 95, label.x=1000, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class)

#Plot coverage controlled background percentages
#Save as 4x16
ggplot(percent_motifs_coverageControl_perSubfamily, aes(Kimura_div, Percent_motifs, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  stat_smooth(method='lm') +
  labs(x='Kimura % divergence', y='Percent elements with motif') +
  stat_cor(label.y = 95, label.x=20, aes(label = ..rr.label..), vjust=-0.5) + 
  stat_cor(label.y = 95, label.x=20, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class, nrow=1)
ggplot(percent_motifs_coverageControl_perSubfamily, aes(Class, Percent_motifs, fill=Class)) +
  geom_violin(bw=10) +
  geom_boxplot(width=0.3) +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  labs(x='', y='Percent elements')
library(dunn.test)
dunn.test(percent_motifs_coverageControl_perSubfamily$Percent_motifs, percent_motifs_coverageControl_perSubfamily$Class, method='bh')
ggplot(percent_motifs_coverageControl_perSubfamily, aes(Length, Percent_motifs, color=Class)) +
  geom_point() +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                          "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,100)) +
  scale_x_log10() +
  stat_smooth(method='lm') +
  labs(x='Consensus length', y='Percent elements with motif') +
  #stat_cor(label.y = 95, label.x=1000, aes(label = ..rr.label..), vjust=-0.5) + 
  #stat_cor(label.y = 95, label.x=1000, aes(label = ..p.label..), vjust=1) +
  facet_wrap(~Class)

