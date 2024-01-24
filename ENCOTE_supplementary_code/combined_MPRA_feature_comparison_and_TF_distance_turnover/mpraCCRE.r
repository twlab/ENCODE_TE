### Monday October 23rd, 2023
# This figure uses the lentivirus MPRA data, cCRE annotations, and TE annotations to measure how different TEs are at driving enhancer activity
# compared to non-TE sequence split by cCRE annotation

# Load in all the necessary libraries
library(ggplot2); library(tidyr)

# Code taken from: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2 
# Creates split violin plots where they are on the same x-axis but each half is a different variable 
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# Now load in the cCRE IDs so you know which ones are TE-derived or not
ccre_ids <- read.table("ccre_IDs.txt"); ccre_te_ids <- read.table("uniq_te_ccre_IDs.txt")
ccre_te_ids <- separate(ccre_te_ids, col = V3, into = c("repClass", "repFam"), sep = "/")
ccre_te_ids$id <- ccre_te_ids$V1; ccre_te_ids <- separate(ccre_te_ids, col = V1, into = c("coords", "ccre"), sep = "\\|"); ccre_te_ids$type <- ifelse(grepl("bound", ccre_te_ids$ccre), "CTCF", ccre_te_ids$ccre)

# The ccre object contains all ccre along with their MPRA value 
ccre <- read.table("k562_cCRE_phastCons_MPRA_bool.txt"); ccre$ID <- ccre$V4; ccre <- separate(ccre, col = V4, into = c("coords", "ccre"), sep = "\\|"); ccre$type <- ifelse(grepl("bound", ccre$ccre), "CTCF", ccre$ccre)
ccre$te <- ifelse(ccre$ID %in% ccre_te_ids$id, "TE", "NonTE")

# Load in the famStats object to find out which TEs are ERV-int
famStats <- read.table("h_famStats.txt"); colnames(famStats)[1] <- "repName"
rownames(famStats) <- famStats$V1; famStats <- separate(famStats, col = V2, into = c("repClass", "repFamily"), sep = "\\/")
internal <- famStats[grepl("-int", rownames(famStats)) | grepl("-I$", rownames(famStats)),"repName"]

# For each of the TE class plots, add in another column with the repClass so you can facet_wrap that. Takes a while...
rownames(ccre_te_ids) <- ccre_te_ids$id; 
ccre_te_ids$repClass <- ifelse(ccre_te_ids$V2 %in% internal, "ERV-int", ifelse(ccre_te_ids$repClass %in% c("LINE", "SINE", "LTR", "DNA"), ccre_te_ids$repClass, "Other"))
ccre$repClass <- ccre_te_ids[ccre$ID, "repClass"]; ccre$repClass <- gsub("\\?", "", ccre$repClass)

# Load in the non-cCRE and non-TE MPRA data
nonAll <- read.table("filt_noncCRE_nonTE_MPRA_bool.txt")

# The TEs which have MPRA data, but are not overlapped with cCREs
nonCCRE_te <- read.table("filt_noncCRE_TE_phastCons_MPRA_bool.txt"); nonCCRE_te <- separate(nonCCRE_te, col = V6, into = c("repClass", "repFam"), sep = "/"); nonCCRE_te$repClass <- gsub("\\?", "", nonCCRE_te$repClass)
nonCCRE_te$repClass <- ifelse(nonCCRE_te$V5 %in% internal, "ERV-int", ifelse(nonCCRE_te$repClass %in% c("LINE", "SINE", "LTR", "DNA"), nonCCRE_te$repClass, "Other") ) 

# Remove the "Low-DNase" category, and create a combined object with all of the different cCRE types for plotting
tmp <- subset(ccre, ccre != "Low-DNase")
full_plot <- data.frame(ccre = c(tmp$type, rep("None", nrow(nonCCRE_te)), rep("None", nrow(nonAll))), mpra = c(tmp$V6, nonCCRE_te$V11, nonAll$V7), bool = c(tmp$V7, nonCCRE_te$V12, rep(0, nrow(nonAll))))
full_plot$te <- c(tmp$te, rep("TE", nrow(nonCCRE_te)), rep("NonTE", nrow(nonAll))) 
full_plot$ccre <- factor(full_plot$ccre, levels = c("None", "DNase-only", "CTCF", "pELS", "dELS", "DNase-H3K4me3", "PLS"))

# Generate objects that contain the total number of annotated elements in each category for geom_text
fL <- as.data.frame(table(subset(full_plot, mpra != 0 & bool != 1)$ccre, subset(full_plot, mpra != 0 & bool != 1)$te)); colnames(fL) <- c("ccre", "te", "Val")
fL_te <- subset(fL, te == "TE"); fL_non <- subset(fL, te == "NonTE")

# Generate the statistical calculations between the TE and nonTE categories
sigs <- data.frame("ccre" = factor(levels(full_plot$ccre), levels = levels(full_plot$ccre)), "mpra" = 0, "te" = "TE", "pval" = 0)
for (f in 1:nrow(sigs)) {
	tmp <- subset(full_plot, mpra != 0 & bool != 1 & ccre == sigs[f,"ccre"])
	sigs[f,"pval"] <- wilcox.test(subset(tmp, te == "TE")$mpra, subset(tmp, te == "NonTE")$mpra, alternative = "greater")$p.value
	sigs[f,"mpra"] <- max(tmp$mpra) + .1
}
# Run p.adjust to correct for multiple tests
sigs$pvalAdj <- p.adjust(sigs$pval); rm(tmp)

# For the significance, you wanted four levels = NS for x > .05, * for .05 > x > .005, ** for .005 > x > .0005, and *** for .0005 > x
sigs$sig <- ifelse(sigs$pvalAdj > .05, "NS", ifelse(sigs$pvalAdj > .005, "*", ifelse(sigs$pvalAdj > .005, "**", "***" ) ) )

# Figure 5a - Split Violin of MPRA data for TEs over cCRE
pdf(paste0(format(Sys.Date(),format="%m%d%y"), "_Fig5b_splitViolin.pdf"))
ggplot(subset(full_plot, mpra != 0 & bool != 1), aes(ccre, mpra, fill = te)) + geom_split_violin() + theme_classic() + labs(x = "cCRE Type", y = "Average Log2FC from MPRA", fill = "Origin") + geom_text(data = fL_non, aes(y = -1.25, label = Val), size = 3, angle = 45, vjust = -2, color = "#F8766D") + geom_text(data = fL_te, aes(y = -1.25, label = Val), size = 3, angle = -45, vjust = -2, color = "#00BFC4") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = .25, position = position_dodge(width = .25), color = "white") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_text(data = sigs, aes(x = ccre, y = mpra, label = sig))
dev.off()



