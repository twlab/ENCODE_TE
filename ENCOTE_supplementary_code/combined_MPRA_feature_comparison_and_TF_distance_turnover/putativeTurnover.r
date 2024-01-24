### Monday October 23rd, 2023
# Load in the files "v2_topMouse_putativeTurnover_phastCons_filter_median" and "v2_topHuman_putativeTurnover_phastCons_filter_median" to 
# make figure 4d that visualizes putative turnover events derived from TEs 

# Load in all the necessary libraries
library(dplyr); library(ggplot2); library(RColorBrewer); library(tidyr); library(scales); library(stringr)

# Now to merge all the colors together into a single object
cols <- c("DNA"="#F8766D", "LINE"="#A3A500", "LTR"="#00BF7D", "ERV-int"="#E76BF3", "SINE"="#00B0F6")

# Old colors...
# cols <- c("LINE" = "#1B9E77", "SINE" = "#D95F02", "LTR" = "#7570B3", "ERV-int" = "#66A61E", "DNA" = "#E7298A")

# Load in the two putative turnover files
hwin5 <- read.table("v2_topHuman_putativeTurnover_phastCons_filter_median", header = T); hwin5$species <- "human"
mwin5 <- read.table("v2_topMouse_putativeTurnover_phastCons_filter_median", header = T); mwin5$species <- "mouse"

# Combine the two objects together and filter
win5 <- rbind(hwin5, mwin5); win5 <- subset( win5, Ancestral_call == "Syntenic" ) 

# Make a row to indicate TE or nonTE
win5$te <- ifelse(win5$TE_ID == ".", 0, 1)

# Generate the base file to plot the TE fraction of putative turnover events
tmp <- subset(win5, te == 1) %>% separate(col = TE_ID, into = c(NA, "repName", "tmp"), sep = "\\|") %>% separate(col = "tmp", into = c("repComb", NA), sep = ",") %>% separate(col = "repComb", into = c("repClass", "repFamily"), sep = "/")

# Change over all the repClass values to be the main 5
internal <- unique(tmp[grepl("-int", tmp$repName) | grepl("-I$", tmp$repName),"repName"])
tmp$repClass <- gsub("\\?", "", tmp$repClass); tmp$repFamily <- gsub("\\?", "", tmp$repFamily); tmp$repClass <- ifelse(tmp$repName %in% internal, "ERV-int", tmp$repClass)
tmp <- subset(tmp, repClass %in% c("LINE", "SINE", "LTR", "ERV-int", "DNA"))

# With the pre-filtering done, time to make the plotting table and make them plots!
put_te_to_plot <- as.data.frame(table(tmp$factor, tmp$species, tmp$repClass)) 
put_te_to_plot <- subset(put_te_to_plot, Freq != 0)
put_te_to_plot$Var3 <- factor(put_te_to_plot$Var3, levels = c("DNA", "LINE", "LTR", "ERV-int", "SINE"))

tot_to <- as.data.frame(table(win5$factor, win5$species)); rownames(tot_to) <- paste0(tot_to$Var1, "|", tot_to$Var2)
put_te_to_plot$Frac <- 0
for (f in 1:nrow(put_te_to_plot)) {
	put_te_to_plot[f,"Frac"] <- round(put_te_to_plot[f,"Freq"] / tot_to[paste0(put_te_to_plot[f,"Var1"], "|", put_te_to_plot[f,"Var2"]),"Freq"], 4)
}
tmp <- subset(put_te_to_plot, Var2 == "human") %>% group_by(Var1) %>% summarize(val = sum(Frac)) %>% arrange(desc(val)) %>% select(Var1) %>% as.data.frame()
put_te_to_plot$Var1 <- factor(put_te_to_plot$Var1, levels = tmp$Var1 )

pdf(paste0(format(Sys.Date(),format="%m%d%y"), "_Fig4d_verticalBarplot.pdf"))
ggplot() + geom_bar(data = subset(put_te_to_plot, Var2 == "human"), aes(fill = Var3, y = Frac, x = Var1), position = "stack", stat = "identity") + geom_bar(data = subset(put_te_to_plot, Var2 == "mouse" ), aes(fill = Var3, y = Frac * -1, x = Var1), position = "stack", stat = "identity") + scale_fill_manual(values = cols, name = "TE Class") + scale_y_continuous(breaks = seq(-.6,.6,.15), labels = c( "60", "45", "30", "15", "0", "15", "30", "45", "60" ), limits = c(-.61,.61) ) + geom_hline(yintercept = 0) + theme_classic() + theme(plot.title = element_text(size = 16, face = "bold.italic"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold")) + labs(x = "", y = "% Putative Turnover Events") + coord_fixed(ratio = 15)
dev.off()
