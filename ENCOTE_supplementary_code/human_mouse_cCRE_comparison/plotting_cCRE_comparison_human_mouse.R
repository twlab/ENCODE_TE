#R 3.5.1
library(ggplot2)
library(XNomial)
library(plyr)
library(ggsci)
source('reference_files/ggplotting_source', encoding = 'UTF-8')

##Orthologous TEs vs. non-TE background
#Human-mouse comparison
#Load data
final_summary_human <- read.table("forR_final_summary_human", header=TRUE, quote="\"")
final_summary_human$Comparison <- factor(final_summary_human$Comparison, levels=c("human_only", "different", "same"))
final_summary_human$cCRE_type <- factor(final_summary_human$cCRE_type, 
                                       levels=c('dELS', 'pELS', 'PLS', 'DNase-H3K4me3', 'CTCF-only'))

#Plot
ggplot(final_summary_human, aes(Annotation, Number, fill=cCRE_type)) +
  geom_bar(position='fill', stat='identity') +
  facet_grid(~Comparison) +
  labs(x='', y='Proportion')
ggplot(final_summary_human, aes(Annotation, Perc_comparison, fill=Comparison)) +
  geom_col() +
  geom_text(aes(label=Number, y=Position)) +
  facet_grid(~cCRE_type) +
  scale_fill_npg() +
  labs(x='', y='Percentage')
#Below plot not very useful
same_only <- subset(final_summary_human, Comparison=='same')
ggplot(same_only, aes(Comparison, Number, fill=Annotation)) +
  geom_bar(position='fill', stat='identity') +
  labs(x='', y='Proportion')

#Mouse-human comparison
#Load data
final_summary_mouse <- read.table("reciprocal_mouse/forR_final_summary_mouse", header=TRUE, quote="\"")
final_summary_mouse$Comparison <- factor(final_summary_mouse$Comparison, levels=c("mouse_only", "different", "same"))
final_summary_mouse$cCRE_type <- factor(final_summary_mouse$cCRE_type, 
                                       levels=c('dELS', 'pELS', 'PLS', 'DNase-H3K4me3', 'CTCF-only'))

ggplot(final_summary_mouse, aes(Annotation, Number, fill=cCRE_type)) +
  geom_bar(position='fill', stat='identity') +
  facet_grid(~Comparison) +
  labs(x='', y='Proportion')
ggplot(final_summary_mouse, aes(Annotation, Perc_comparison, fill=Comparison)) +
  geom_col() +
  geom_text(aes(label=Number, y=Position)) +
  facet_grid(~cCRE_type) +
  scale_fill_npg() +
  labs(x='', y='Percentage')

##Novel cCREs from TEs vs. non-TE background
#Human novel cCREs
#Load data
summary_novel_human <- read.delim("forR_summary_novel_human", comment.char="#")

#Plot
ggplot(summary_novel_human, aes(Annotation, Perc_TE, fill=cCRE_type)) +
  geom_bar(position='fill', stat='identity') +
  labs(x='', y='Proportion')
ggplot(summary_novel_human, aes(cCRE_type, Perc_cCRE, fill=Annotation)) +
  geom_bar(position='fill', stat='identity') +
  labs(x='', y='Proportion')

#Mouse novel cCREs
#Load data
summary_novel_mouse <- read.delim("reciprocal_mouse/forR_summary_novel_mouse", comment.char="#")

#Plot
ggplot(summary_novel_mouse, aes(Annotation, Perc_TE, fill=cCRE_type)) +
  geom_bar(position='fill', stat='identity') +
  labs(x='', y='Proportion')
ggplot(summary_novel_mouse, aes(cCRE_type, Perc_cCRE, fill=Annotation)) +
  geom_bar(position='fill', stat='identity') +
  labs(x='', y='Proportion')

#Combine human and mouse novel cCREs
summary_novel_human$Species <- "Human"
summary_novel_mouse$Species <- "Mouse"
novel_human_mouse_ccres_combined <- rbind(summary_novel_human, summary_novel_mouse)
novel_human_mouse_ccres_combined <- subset(novel_human_mouse_ccres_combined,
                                           Annotation != 'non-TE')
novel_human_mouse_ccres_combined$cCRE_type <- factor(novel_human_mouse_ccres_combined$cCRE_type,
                                                     levels=c('dELS', 'pELS', 'PLS', 'DNase-H3K4me3', 'CTCF-only'))
ggplot(novel_human_mouse_ccres_combined, aes(cCRE_type, Perc_cCRE)) +
  stat_summary(fun.data = mean_sdl, geom='bar') +
  geom_jitter(aes(color=Species), width=0, height=0) +
  labs(x='', y='Percentage')

#Multinomial exact test for TE vs non-TE (Human)
#PLS
human_PLS_TE_num <- c(unlist(subset(final_summary_human, Annotation=='TE' & cCRE_type=='PLS', select=Number)))
human_PLS_nonTE_num <- c(unlist(subset(final_summary_human, Annotation=='non-TE' & cCRE_type=='PLS', select=Number)))
xmulti(human_PLS_TE_num, human_PLS_nonTE_num, detail=3, histobins = T) #LLR p-value 3.922e-26
#pELS
human_pELS_TE_num <- c(unlist(subset(final_summary_human, Annotation=='TE' & cCRE_type=='pELS', select=Number)))
human_pELS_nonTE_num <- c(unlist(subset(final_summary_human, Annotation=='non-TE' & cCRE_type=='pELS', select=Number)))
xmonte(human_pELS_TE_num, human_pELS_nonTE_num, ntrials=1000000, detail=3, histobins = T) #LLR p-value 0
#dELS
human_dELS_TE_num <- c(unlist(subset(final_summary_human, Annotation=='TE' & cCRE_type=='dELS', select=Number)))
human_dELS_nonTE_num <- c(unlist(subset(final_summary_human, Annotation=='non-TE' & cCRE_type=='dELS', select=Number)))
xmonte(human_dELS_TE_num, human_dELS_nonTE_num, ntrials=1000000, detail=3, histobins = T) #LLR p-value 0
#DNase-H3K4me3
human_H3K4me3_TE_num <- c(unlist(subset(final_summary_human, Annotation=='TE' & cCRE_type=='DNase-H3K4me3', select=Number)))
human_H3K4me3_nonTE_num <- c(unlist(subset(final_summary_human, Annotation=='non-TE' & cCRE_type=='DNase-H3K4me3', select=Number)))
xmulti(human_H3K4me3_TE_num, human_H3K4me3_nonTE_num, detail=3, histobins = T) #LLR p-value 3.975e-16
#CTCF-only
human_CTCF_TE_num <- c(unlist(subset(final_summary_human, Annotation=='TE' & cCRE_type=='CTCF-only', select=Number)))
human_CTCF_nonTE_num <- c(unlist(subset(final_summary_human, Annotation=='non-TE' & cCRE_type=='CTCF-only', select=Number)))
xmulti(human_CTCF_TE_num, human_CTCF_nonTE_num, detail=3, histobins = T) #LLR p-value 7.233e-28

#Compare exact multinomial results to G-test (Human)
library(AMR)
g.test(PLS_TE_num, p=PLS_nonTE_num, rescale.p = TRUE) #p-value < 2.2e-16
g.test(pELS_TE_num, p=pELS_nonTE_num, rescale.p = TRUE) #p-value < 2.2e-16
g.test(dELS_TE_num, p=dELS_nonTE_num, rescale.p = TRUE) #p-value < 2.2e-16
g.test(H3K4me3_TE_num, p=H3K4me3_nonTE_num, rescale.p = TRUE) #p-value = NA (probably because of the 0 in observed)
g.test(CTCF_TE_num, p=CTCF_nonTE_num, rescale.p = TRUE) #p-value < 2.2e-16

#Multinomial exact test for TE vs non-TE (Mouse)
#PLS
mouse_PLS_TE_num <- c(unlist(subset(final_summary_mouse, Annotation=='TE' & cCRE_type=='PLS', select=Number)))
mouse_PLS_nonTE_num <- c(unlist(subset(final_summary_mouse, Annotation=='non-TE' & cCRE_type=='PLS', select=Number)))
xmulti(mouse_PLS_TE_num, mouse_PLS_nonTE_num, detail=3, histobins = T) #LLR p-value 2.008e-10
#pELS
mouse_pELS_TE_num <- c(unlist(subset(final_summary_mouse, Annotation=='TE' & cCRE_type=='pELS', select=Number)))
mouse_pELS_nonTE_num <- c(unlist(subset(final_summary_mouse, Annotation=='non-TE' & cCRE_type=='pELS', select=Number)))
xmonte(mouse_pELS_TE_num, mouse_pELS_nonTE_num, ntrials=1000000, detail=3, histobins = T) #LLR p-value 0
#dELS
mouse_dELS_TE_num <- c(unlist(subset(final_summary_mouse, Annotation=='TE' & cCRE_type=='dELS', select=Number)))
mouse_dELS_nonTE_num <- c(unlist(subset(final_summary_mouse, Annotation=='non-TE' & cCRE_type=='dELS', select=Number)))
xmonte(mouse_dELS_TE_num, mouse_dELS_nonTE_num, ntrials=1000000, detail=3, histobins = T) #LLR p-value 0
#DNase-H3K4me3
mouse_H3K4me3_TE_num <- c(unlist(subset(final_summary_mouse, Annotation=='TE' & cCRE_type=='DNase-H3K4me3', select=Number)))
mouse_H3K4me3_nonTE_num <- c(unlist(subset(final_summary_mouse, Annotation=='non-TE' & cCRE_type=='DNase-H3K4me3', select=Number)))
xmulti(mouse_H3K4me3_TE_num, mouse_H3K4me3_nonTE_num, detail=3, histobins = T) #LLR p-value 0.001799
#CTCF-only
mouse_CTCF_TE_num <- c(unlist(subset(final_summary_mouse, Annotation=='TE' & cCRE_type=='CTCF-only', select=Number)))
mouse_CTCF_nonTE_num <- c(unlist(subset(final_summary_mouse, Annotation=='non-TE' & cCRE_type=='CTCF-only', select=Number)))
xmulti(mouse_CTCF_TE_num, mouse_CTCF_nonTE_num, detail=3, histobins = T) #LLR p-value 2.645e-20

################### Novel and conserved TE contributions ####################
#Load data
#File put together manually from final_summary_human and final_summary_mouse files
novel_and_conserved_TE_contributions <- read.delim("novel_and_conserved_TE_contributions")
novel_and_conserved_TE_contributions_TEonly <- subset(novel_and_conserved_TE_contributions,
                                                      Annotation == 'TE')
novel_and_conserved_TE_contributions_TEonly$cCRE_type <- factor(novel_and_conserved_TE_contributions_TEonly$cCRE_type,
                                                     levels=c('dELS', 'pELS', 'PLS', 'DNase-H3K4me3', 'CTCF-only'))
#Plot
ggplot(novel_and_conserved_TE_contributions_TEonly, aes(Function, Perc_category)) +
  stat_summary(fun.data = mean_sdl, geom='bar') +
  geom_jitter(aes(color=Species), width=0, height=0) +
  labs(x='', y='Percentage') +
  facet_grid(cols=vars(cCRE_type))

