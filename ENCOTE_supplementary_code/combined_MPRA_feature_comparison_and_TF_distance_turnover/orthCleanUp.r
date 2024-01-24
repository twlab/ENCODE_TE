### Monday March 13th, 2022
# This code will be invoked with 'Rscript' and will load in a list of peaks that contain a column of IDs (some with multiple separated by ',') and another
# column of integers (also ',' separated). Turn those separated values into a list and only keep the one associated with the larger value (or the 1st if 
# there's a tie). 
#
# Expected script call (post running module load R/4.1.3-r9):
# Rscript --vanilla orthCleanUp.r [peak file] [factor output file] [second run?]

# Load in the cool package stringr
library(stringr)

args = commandArgs(trailingOnly = T)

# Load in the file
b <- read.table(args[1])
t <- subset(b, grepl(",", V4))

i <- str_split(t$V4, ",")
v <- str_split(t$V5, ",")

# This pulls out (as an integer array) which index is the largest (first if ties) across each list
h <- sapply(v, function(x) which.max(x) )

# Use a good ol' fashioned loop to iterate through the ids and replace the value in the original file
for (f in 1:length(h)) { b[rownames(t)[f], "V4"] <- i[[f]][h[f]]; b[rownames(t)[f], "V5"] <- v[[f]][h[f]] }

# With everything done, time to output the finished product, yay
if (args[3] == "F") { write.table(b, file = args[2], quote = F, row.names = F, col.names = F, sep = "\t") } else { write.table(b[,c(1:4)], file = args[2], quote = F, row.names = F, col.names = F, sep = "\t") }
