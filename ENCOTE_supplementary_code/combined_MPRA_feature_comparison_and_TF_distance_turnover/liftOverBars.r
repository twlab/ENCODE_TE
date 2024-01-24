### Monday October 23rd, 2023
# This code will load in the liftOver files and combine them together to make a vertical barplot for figure 4c

# Load in required libraries
library(dplyr); library(ggplot2); library(RColorBrewer); library(tidyr); library(scales); library(stringr); require(doParallel)

# Load in the paths/files you need for the analysis
samps <- read.table("fullFactor")$V1

# Was originally pulling in the "_full" files, which shouldn't be much different...?

# Load in all the files and add additional columns to specify factor and species
results_orth <- mclapply(1:length(samps), function(i) {
	tmpM <- read.table(paste0("1_mouseTFs/", samps[i], "_full_5kb_merge"), header = T); tmpM$sample <- samps[i]; tmpM$species <- "mouse"; rownames(tmpM) <- tmpM$newSpeciesID
	tmpH <- read.table(paste0("2_humanTFs/", samps[i], "_full_5kb_merge"), header = T); tmpH$sample <- samps[i]; tmpH$species <- "human"; rownames(tmpH) <- tmpH$newSpeciesID
	rbind(subset(tmpH, grepl("\\|H", newSpeciesID)), subset(tmpM, grepl("\\|M", newSpeciesID)))
}, mc.cores = 5)

# Combine all of the values together into a single data.frame
orth <- bind_rows(results_orth)

# Remove a lot of the non-TE stuff since you want to focus on TEs and their contribution to orthology and turnover events! 
out_plot <- orth; out_plot <- out_plot[rowSums(out_plot == -2) == 0,] 

# In order to keep things identical to how it was done before, change mate_same_te to 0 if teID does not equal "." and isn't already marked as "1"
out_plot$mate_same_te <- ifelse(out_plot$teID != "." & out_plot$mate_same_te != 1, 0, out_plot$mate_same_te)
out_tab_plot <- as.data.frame(table(out_plot$sample, out_plot$orth, out_plot$mate_act, out_plot$mate_same_te, out_plot$species))
out_tab_plot$Qual <- paste(out_tab_plot$Var2, out_tab_plot$Var3, out_tab_plot$Var4, sep = ":"); out_tab_plot <- subset(out_tab_plot, Freq != 0)
out_tab_plot$Qual <- factor(out_tab_plot$Qual, levels = c("0:0:-1", "0:0:0", "1:0:-1", "1:1:-1", "1:0:0", "1:1:0", "1:1:1"))

# Now add total fraction to each class to make things easier
tot <- as.data.frame(table(orth$sample, orth$species)); rownames(tot) <- paste0(tot$Var1, "|", tot$Var2)
out_tab_plot$Frac <- 0
for (f in 1:nrow(out_tab_plot)) {
	out_tab_plot[f,"Frac"] <- round(out_tab_plot[f,"Freq"] / tot[paste0(out_tab_plot[f,"Var1"], "|", out_tab_plot[f,"Var5"]),"Freq"], 4)
}

# Combine all of the various TE information together to simplify things
out_tab_plot$Qual <- as.character(out_tab_plot$Qual); out_tab_plot$Qual <- ifelse( out_tab_plot$Qual == "1:0:1", "1:0:0", ifelse( out_tab_plot$Qual == "1:1:1", "1:1:0", out_tab_plot$Qual ) )
tmp <- out_tab_plot[out_tab_plot$Var5 != "mouse" & out_tab_plot$Var4 != -1,] %>% group_by(Var1) %>% summarize(Frac = sum(Frac)) %>% arrange(desc(Frac)) %>% select(Var1) %>% as.data.frame()
out_tab_plot$Var1 <- factor(out_tab_plot$Var1, levels = tmp$Var1)

# Generate figure 4c
valsNew <- c("0:0:0" = "Species-Specific", "1:0:0" = "Syntenic but Specific", "1:1:0" = "Syntenic-Shared")
pdf(paste0(format(Sys.Date(),format="%m%d%y"), "_Fig4c_verticalBarplot.pdf"))
ggplot() + geom_bar(data = subset(out_tab_plot, Var5 == "human" & Qual %in% c("0:0:0", "1:0:0", "1:1:0") ), aes(fill = Qual, y = Frac, x = Var1), position = "stack", stat = "identity") + geom_bar(data = subset(out_tab_plot, Var5 == "mouse" & Qual %in% c("0:0:0", "1:0:0", "1:1:0") ), aes(fill = Qual, y = Frac * -1, x = Var1), position = "stack", stat = "identity") + scale_fill_brewer(palette = "Paired", labels = valsNew) + scale_y_continuous(breaks = seq(-.45,.45,.15), labels = c( "45", "30", "15", "0", "15", "30", "45" ), limits = c(-.45,.45) ) + geom_hline(yintercept = 0) + theme_classic() + theme(plot.title = element_text(size = 16, face = "bold.italic"), axis.text = element_text(size = 12), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold")) + labs(x = "", y = "Percentage of TF Binding Sites", fill = "TE Category") + coord_flip()
dev.off()

