#R 3.5.1
library(ggplot2)
source('reference_files/ggplotting_source', encoding = 'UTF-8')

#####TE family level enrichment for cCRE overlapping regions over foreground sequence#####
TEfamily_ccreOnly_overlap_coverage <- read.delim("ccreOnly_overlap_coverage/TEfamily_ccreOnly_overlap_coverage")
TEfamily_ccreOnly_overlap_coverage <- subset(TEfamily_ccreOnly_overlap_coverage,
                                             Family != "DNA/TcMar-Tc1" & Family != "LINE/Dong-R4" & Family != "LINE/Jockey"
                                             & Family != "LINE/Penelope" & Family != "Other")
ggplot(TEfamily_ccreOnly_overlap_coverage, aes(Position, Enrichment, color=Family)) +
  geom_line() +
  geom_hline(yintercept = 1, color='grey', linetype='dashed') +
  #scale_y_continuous(limits=c(0,3)) +
  facet_wrap(~Family) +
  theme(legend.position="none") +
  theme(plot.margin = margin(10, 20, 10, 10)) +
  labs(x='Normalized consensus position', y='Relative enrichment')

#####TE class and family level enrichment for cCRE vs. non-cCRE#####
combined_normEnrichment_TEclass <- read.delim("length_controlled_random_background_coverage_max_v3/combined_normEnrichment_TEclass")
combined_normEnrichment_TEfamily <- read.delim("length_controlled_random_background_coverage_max_v3/combined_normEnrichment_TEfamily")
combined_normEnrichment_TEfamily_wClass <- read.delim("length_controlled_random_background_coverage_max_v3/combined_normEnrichment_TEfamily_wClass")
#Plot class level enrichment
combined_normEnrichment_TEclass <- subset(combined_normEnrichment_TEclass, Class != "Other" & Class != "SVA")
ggplot(combined_normEnrichment_TEclass, aes(Consensus_bin, Enrichment, color=Class)) +
  geom_line() +
  geom_hline(yintercept = 1, color='grey', linetype='dashed') +
  geom_ribbon(aes(x=Consensus_bin, ymin=Lower, ymax=Upper, fill=Class), alpha=0.2) +
  scale_fill_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_color_manual("TE class", values=c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D",
                                         "SINE"="#00B0F6", "ERV-int"="#E76BF3")) +
  scale_y_continuous(limits=c(0,4)) +
  facet_wrap(~Class, nrow=1) +
  theme(legend.position="none") +
  labs(x='Normalized consensus position', y='Relative enrichment')
#Plot family level enrichment
ggplot(combined_normEnrichment_TEfamily, aes(Consensus_bin, Enrichment, color=Family)) +
  geom_line() +
  geom_hline(yintercept = 1, color='grey', linetype='dashed') +
  geom_ribbon(aes(x=Consensus_bin, ymin=Lower, ymax=Upper, fill=Family), alpha=0.2) +
  scale_y_continuous(limits=c(0,4)) +
  facet_wrap(~Family) +
  theme(legend.position="none") +
  theme(plot.margin = margin(10, 20, 10, 10)) +
  labs(x='Normalized consensus position', y='Relative enrichment')
#Plot family level enrichment with additional TE class info
ggplot(combined_normEnrichment_TEfamily_wClass, aes(Consensus_bin, Enrichment, color=Family)) +
  geom_line() +
  geom_hline(yintercept = 1, color='grey', linetype='dashed') +
  geom_ribbon(aes(x=Consensus_bin, ymin=Lower, ymax=Upper, fill=Family), alpha=0.2) +
  scale_y_continuous(limits=c(0,4)) +
  facet_wrap(~Family) +
  theme(legend.position="none") +
  theme(plot.margin = margin(10, 20, 10, 10)) +
  labs(x='Normalized consensus position', y='Relative enrichment')

