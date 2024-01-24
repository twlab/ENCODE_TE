### Sunday October 22nd, 2023
# This code will load in the given peak file for each downloaded hg38 ChIP-seq factor and the hg38 TE annotations 

# Load in required libraries
library(data.table); library(dplyr); library(regioneR); library(stringi); require(doParallel); library(ggplot2); library(tidyr); library(scales)

# Load in all the pre-intercepted files between TEs and peaks  
ids <- list.files("3_distance/", "_IDs$"); samps <- gsub(".{4}$", "", ids)

# Set this so that when you use fread it doesn't print an annoying message
options(datatable.fread.input.cmd.message=FALSE)

# Load in the TE annotations
TEs <- fread("TE_hg38.txt.gz"); colnames(TEs) <- c("chr", "start", "end", "repName", "repComb", "strand", "consensus", "div", "ID")
TEs <- makeGRangesFromDataFrame(subset(TEs, chr != "chrY"), keep.extra.columns=T, starts.in.df.are.0based=F); seqlevels(TEs) <- sort(seqlevels(TEs))

# Code taken from: https://support.bioconductor.org/p/110506/ and is the 3rd version
# This will return an output similar to bedtools closest -io which ignores overlaps for the nearest element
fun2 <- function(x) min(x) == x & !duplicated(x)
dTN <- function(query, subject, fun = fun2) {
	hits <- union(precede(query, subject, select="all"), follow(query, subject, select="all"))
	distance <- distance(query[queryHits(hits)], subject[subjectHits(hits)])
	lst <- splitAsList(distance, queryHits(hits))
	keep <- unsplit(fun(lst), queryHits(hits))
	mcols(hits)$distance <- distance
	hits[keep]
}

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


# Remove all the samples with *no* TE families with at least 10 bound elements and keep the various file lists in the same order
filt <- !(file.size(paste0("3_distance/", samps, "_10fams")) == 0)
samps <- samps[filt]; ids <- ids[filt] 


set.seed(8)
res <- mclapply(1:length(samps), function(i) {

	# Load in the file for each sample
	tmp <- read.table(paste0("3_distance/", ids[i]))

	# Set the column names so everything is consistent rather than V#
	colnames(tmp) <- c("chr", "start", "end", "sumStart", "sumEnd", "peakID", "teChr", "teStart", "teEnd", "repName", "repComb", "strand", "cons", "div", "teID") 
	tmp <- makeGRangesFromDataFrame(subset(tmp, chr != "chrY"), keep.extra.columns=T, starts.in.df.are.0based=T); seqlevels(tmp) <- sort(seqlevels(tmp))
	
	# Make a table of all of the families that have at least 10 elements intersecting the particular factor
	tmp_tbl <- data.frame(table(gsub(".*\\|", "", unique(subset(tmp, teID != ".")$teID))))
	tmp_fams <- subset(tmp_tbl, Freq >= 10); rownames(tmp_fams) <- tmp_fams$Var1

	# Make a data.frame to store all the results into 
	df <- data.frame(samp = samps[i], repName = rownames(tmp_fams)[1]); df$within <- list(1); df$without_mean <- list(1) 
	df$without_rank <- list(1); df$pVal_mean <- 1; df$pVal_rank <- 1; df$pVal_all <- 1; df$numNotBound <- 1; df$notBoundSummary <- list(1)

	# Do not perform the analysis if there are no families with more than 10 elements 
	if (nrow(tmp_fams) >= 1) {
	
		# Make a copy of TEs so that you can append the distances for each factor as needed. If there is an error (such as a factor 
		# only having one peak on a chromosome and it overlaps a TE, making that TE have no distance value) then return a filter value
		tmp_te <- TEs; dist <- as.data.frame(dTN(tmp_te, tmp))$distance
		if (length(dist) != length(tmp_te)) {
			df$repName <- "FAIL"; return(df)
		} else { 
			tmp_te$peakDist <- dist
		}

		# Since there are families with at least 10 elements, perform the full analysis
		for (f in 1:nrow(tmp_fams)) {

			if (f == 1) {		
				tmp_df <- df
			} else {
				tmp_df[f,"samp"]<- samps[i]; tmp_df[f,"repName"]<- rownames(tmp_fams)[f]
			}

			# Make a subset of elements for the specific family that you are currently analyzing and mark them as bound or not
			sub <- subset(tmp_te, repName == rownames(tmp_fams)[f]); sub$bound <- ifelse(sub$ID %in% tmp$teID, T, F)

			# Save the distances of elements bound with the factor summit		
			tmp_df[f,"within"][[1]] <- list(subset(sub, bound == T)$peakDist)

			# Now bootstrap drawing tmp_fams$V2 number of non-overlapped elements 1000 times to get a list of mean random distances
			boot <- subset(sub, bound == F)$peakDist
			tmp_list <- c()
			for (a in 1:1000) { tmp_list <- c(tmp_list, mean(boot[sample(1:length(boot), tmp_fams[f,1], T)])) }
			tmp_df[f,"without_mean"][[1]] <- list(tmp_list)
		
			# Now do the same as above, but instead get a rank list of mean values to recapitulate the expected distribution 	
			tmp_list <- c() 
			for (a in 1:1000) { tmp_list[[a]] <- sort(boot[sample(1:length(boot), tmp_fams[f,1], T)]) }
			tmp_df[f,"without_rank"][[1]] <- list(colMeans(do.call(rbind, tmp_list)))

			# Finally, generate the test comparing all bound elements to all unbound elements within the same family
			tmp_df[f,"pVal_all"] <- tryCatch(wilcox.test(unlist(tmp_df[f,"within"]), boot, exact = F)$p.value, error = function(e) 1)
			tmp_df[f,"numNotBound"] <- length(boot)
			tmp_df[f,"notBoundSummary"][[1]] <- list(as.numeric(summary(log10(boot+1))[c(1:3,5:6)]))

			# Now for p-value testing
			tmp_df[f,"pVal_mean"] <- tryCatch(wilcox.test(unlist(tmp_df[f,"within"]), unlist(tmp_df[f,"without_mean"]), exact = F)$p.value, error = function(e) 1)
			tmp_df[f,"pVal_rank"] <- tryCatch(wilcox.test(unlist(tmp_df[f,"within"]), unlist(tmp_df[f,"without_rank"]), exact = F)$p.value, error = function(e) 1)
		}
	}
	# Output the complete data.frame
	return(tmp_df)

}, mc.cores = 60)

# Bind all the rows of the above object to easily parse through the output. Also remove instances where there were not enough peaks
bind <- bind_rows(res); bind <- subset(bind, repName != "FAIL")

# Since you love lists, turn bind back into a list along each sample for the below code
nres <- list(); uniq <- unique(bind$samp)
for (f in 1:length(uniq)) {
	nres[[f]] <- subset(bind, samp == uniq[f])
}

# Make another class for the ERV-int as before
famStats <- read.table("h_famStats.txt"); colnames(famStats)[1] <- "repName"
rownames(famStats) <- famStats$V1; famStats <- separate(famStats, col = V2, into = c("repClass", "repFamily"), sep = "\\/")
internal <- famStats[grepl("-int", rownames(famStats)) | grepl("-I$", rownames(famStats)),"repName"]
famStats$repClass <- gsub("\\?", "", famStats$repClass); famStats$repFamily <- gsub("\\?", "", famStats$repFamily)
famStats$repClass <- ifelse(famStats$repName %in% internal, "ERV-int", famStats$repClass)
rownames(famStats) <- famStats$repName

# To make the above split_violin plot, you'll need to combine all of the data in 'res' that includes the ranked "without" distances and the 
# "within" bound distances  
resp <- famStats; resp$without_rank <- resp$within <- rep(list(list()), nrow(resp))
for (f in 1:length(nres)) {
	for (g in 1:nrow(nres[[f]])) {
		cur <- nres[[f]][g,"repName"]
		resp[cur, "within"][[1]] <- list(c(resp[cur,"within"][[1]], nres[[f]][g,"within"]))
		resp[cur, "without_rank"][[1]] <- list(c(resp[cur,"without_rank"][[1]], nres[[f]][g,"without_rank"]))
	}
}

# Now that you've combined all the results of res, split it by repClass
sine <- subset(resp, repClass == "SINE"); sine$wiMed <- lapply(sine$within, function(x) {lapply(x, function(y) {median(y)} )} ); sine$woMed <- lapply(sine$without_rank, function(x) {lapply(x, function(y) {median(y)} )} )
line <- subset(resp, repClass == "LINE"); line$wiMed <- lapply(line$within, function(x) {lapply(x, function(y) {median(y)} )} ); line$woMed <- lapply(line$without_rank, function(x) {lapply(x, function(y) {median(y)} )} )
ltr <- subset(resp, repClass == "LTR" | repClass == "ERV-int"); ltr$wiMed <- lapply(ltr$within, function(x) {lapply(x, function(y) {median(y)} )} ); ltr$woMed <- lapply(ltr$without_rank, function(x) {lapply(x, function(y) {median(y)} )} )
dna <- subset(resp, repClass == "DNA"); dna$wiMed <- lapply(dna$within, function(x) {lapply(x, function(y) {median(y)} )} ); dna$woMed <- lapply(dna$without_rank, function(x) {lapply(x, function(y) {median(y)} )} )
comb <- rbind(sine, line, ltr, dna)

# Now turn each median value within all of the above lists into a row in a data.frame with a column for "within" and "without" to be reshaped for ggplot
vals <- data.frame()
for (f in c("LINE", "SINE", "DNA", "LTR")) {
	sub <- subset(comb, repClass == f)
	vals <- rbind(vals, data.frame("samp" = rep( f, length(unlist(sub$wiMed)) ), "win" = unlist( sub$wiMed ), "wout_rank" = unlist( sub$woMed ) ) )
}

vals <- reshape2::melt(vals, value.name = "med_nearest_peak_dist"); vals$variable <- ifelse( vals$variable == "win", "bound", "nonBound" )
vals$samp <- factor(vals$samp, levels = c("LINE", "SINE", "LTR", "DNA") )

# Generate the wilcoxon significance vals for each comparison 
sigs <- data.frame("ccre" = factor(levels(vals$samp), levels = levels(vals$samp)), "med_nearest_peak_dist" = 0, "variable" = "bound", "pval" = 0)
for (f in 1:nrow(sigs)) {
	tmp <- subset(vals, samp == sigs[f,1])
	sigs[f,"pval"] <- wilcox.test(subset(tmp, variable == "bound")$med_nearest_peak_dist, subset(tmp, variable == "nonBound")$med_nearest_peak_dist, alternative = "less")$p.value
	sigs[f,"med_nearest_peak_dist"] <- max(tmp$med_nearest_peak_dist) + (.1 * max(tmp$med_nearest_peak_dist))
}
# Now run p.adjust to correct for multiple valss
sigs$pvalAdj <- p.adjust(sigs$pval)

# For the significance, you wanted four levels = NS for x > .05, * for .05 > x > .005, ** for .005 > x > .0005, and *** for .0005 > x
sigs$sig <- ifelse(sigs$pvalAdj > .05, "NS", ifelse(sigs$pvalAdj > .005, "*", ifelse(sigs$pvalAdj > .005, "**", "***" ) ) )

# Finally, make the split violin plot for the median ranked distance
pdf(paste0(format(Sys.Date(),format="%m%d%y"), "_Fig4b_splitViolin.pdf"))
ggplot(data = vals, aes(x = samp, y = med_nearest_peak_dist, fill = variable)) + geom_split_violin() + theme_classic() + theme(axis.text.y = element_text(size = 16), axis.title = element_text(size = 18, face = "bold.italic"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12)) + labs(x = "", y = "Median Distance to Nearest Peak") + scale_fill_discrete(name = "", labels = c("TF bound TEs", "Not bound TEs")) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + annotation_logticks(side = "l") + geom_text(data = sigs, aes(x = ccre, y = med_nearest_peak_dist, label = sig))
dev.off()

