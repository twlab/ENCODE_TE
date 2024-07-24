### Monday October 23rd, 2023
# This figure uses the lentivirus MPRA data, TE annotations, K562 TF ChIP-seq, and PhastCons to measure how human lineage intervals (both TE
# and nonTE) overlap with genomic indicators of cis-regulatory activity

# Load in all the necessary libraries
library(data.table); library("ggVennDiagram"); library(ggplot2)

# Read in the ChIP-seq intersection IDs
chip_te <- subset(read.table("combined_ccre_chip"), V5 == "TE")$V4
chip_nonte <- subset(read.table("combined_ccre_chip"), V5 == "nonTE")$V4

# Read in the ATAC-seq intersection IDs
atac_te <- read.table("ccre_te_atac")$V1
atac_nonte <- read.table("ccre_nonte_atac")$V1

# Read in the MPRA intersection IDs
mpra_te <- read.table("ccre_te_mpra")$V1
mpra_nonte <- read.table("ccre_nonte_mpra")$V1

# Read in the phastCons intersection IDs
bins <- fread("bins_noNa_phastCons30.txt")
phast_te <- subset(bins, V6 >= quantile(bins$V6, .9) & V5 == "TE")$V4
phast_nonte <- subset(bins, V6 >= quantile(bins$V6, .9) & V5 == "nonTE")$V4

# Create the venn diagram objects
# TE
i <- list("TF ChIP-seq" = chip_te, "ATAC-seq" = atac_te, "30-way phastCons" = phast_te, "MPRA" = mpra_te)
te <- Venn(i); te <- process_data(te)

# nonTE
i <- list("TF ChIP-seq" = chip_nonte, "ATAC-seq" = atac_nonte, "30-way phastCons" = phast_nonte, "MPRA" = mpra_nonte)
nonte <- Venn(i); nonte <- process_data(nonte)

# Calculate chi-square test results for TE vs nonTE proportions
te_v <- venn_region(te)$count; nonte_v <- venn_region(nonte)$count
res <- chisq.test(te_v, nonte_v)$p.value; print(paste0("Chi-Square P-value: ", res))	

pdf(paste0(format(Sys.Date(),format="%m%d%y"), "_Fig5c_vennDiagram.pdf"))
ggplot() + geom_sf(data = venn_setedge(te), linewidth = 2.25, aes(color = name), show.legend = F) + geom_sf_label(data = (venn_region(te) %>% dplyr::mutate(percent = paste(round(.data$count * 100/sum(.data$count), 2))) ), aes(label = paste0(percent, "%\n", formatC(count, format="d", big.mark=","))), fontface = "bold") + geom_sf_text(data = venn_setlabel(te), aes(label = name)) + theme_void() + ggtitle("TE Intervals")
ggplot() + geom_sf(data = venn_setedge(nonte), linewidth = 2.25, aes(color = name), show.legend = F) + geom_sf_label(data = (venn_region(nonte) %>% dplyr::mutate(percent = paste(round(.data$count * 100/sum(.data$count), 2))) ), aes(label = paste0(percent, "%\n", formatC(count, format="d", big.mark=","))), fontface = "bold") + geom_sf_text(data = venn_setlabel(nonte), aes(label = name)) + theme_void()  + ggtitle("nonTE Intervals")
dev.off()
