### Sunday March 13th, 2022
# This code will load in the given file(s) for each of the factors manually so that the various objects needed for plotting can be made together
# and combined all within R. 

# Load in required libraries
library(tidyr); library(dplyr); library(data.table); library(stringr); require(doParallel)

# Read in which factors you'll iterate over
samps <- read.table("fullFactor")$V1

results <- mclapply(1:length(samps), function(i) {

	hp <- fread(paste0("2_humanTFs/", samps[i])); mp <- fread(paste0("1_mouseTFs/", samps[i]))
	hlo <- fread(paste0("2_humanTFs/mm10.hg38.", samps[i], ".int2")); mlo <- fread(paste0("1_mouseTFs/hg38.mm10.", samps[i], ".int2"))
	hwin <- read.table(paste0("2_humanTFs/", samps[i], "_5kb")); mwin <- read.table(paste0("1_mouseTFs/", samps[i], "_5kb")) 

	# Set the rownames as needed for ease of data access
	rownames(hp) <- hp$V5; rownames(mp) <- mp$V5

	# Find the IDs which occur multiple times in the liftOver file so you can deal with them separately later. It should be noted that the 
	# hdupIDs are actually IDs from mouse intervals and vice versa for mdupIDs
	hdupIDs <- as.character(subset(arrange(as.data.frame(table(hlo$V5)), desc(Freq)), Freq > 1)$Var1)
	mdupIDs <- as.character(subset(arrange(as.data.frame(table(mlo$V5)), desc(Freq)), Freq > 1)$Var1)

	# Initiallize hout/mout further down for ease of troubleshooting/testing. These files only contain peaks that liftedOver at least once!
	hout <- mout <- data.table(matrix(0, nrow = 0, ncol = 0))
	hout <- data.frame(speciesID = unique(hlo$V5)); mout <- data.frame(speciesID = unique(mlo$V5))
	rownames(hout) <- hout$speciesID; rownames(mout) <- mout$speciesID

	# Set the columns for orthology and mate_act to -1 and parse through in the loop below
	hout$mate_same_te <- hout$mate_act <- hout$orth <- mout$mate_same_te <- mout$mate_act <- mout$orth <- -1

	# Iterate through hlo using the unique IDs as seen in hout and determine whether the particular ID is orthologous or not, as well as 
	# whether the orthologous TE-derived peak is derived from the same TE in both human and mouse or a different one (using repFamily!)
	for (f in 1:nrow(hout)) {
		hich <- which(hlo$V5 == hout[f,1]); mich <- which(mlo$V7 == hout[f,1])
		hsub <- hlo[hich,]; msub <- mlo[mich,]

		### ORTHOLOGY
		# If at least 1 hsub instance reciprocally maps back to itself with the current hout value, orthologous. Otherwise, not orthologous 
		if ( sum( hsub$V5 == hsub$V9 ) > 0 ) { hout[f,"orth"] <- 1 
			} else { hout[f,"orth"] <- hout[f,"mate_act"] <- 0; if ( grepl("\\|TE", hout[f,"speciesID"]) ) { hout[f,"mate_same_te"] <- 0 }  }

		### MATE ACTIVITY
		# If not orthologous, 0 (change nothing as done above). If orthologous, check for at least one peak ID in hsub or no mich length for 0/1, respectively
		if ( hout[f,"orth"] == 1 && ( sum( hsub$V7 != "." ) > 0 || length(mich) > 0 ) ) { hout[f,"mate_act"] <- 1 } else { hout[f,"mate_act"] <- 0 }

		### SAME TE
		# If not orthologous AND TE, 0 (no change). If not TE, -1 (no change). Otherwise, check if the family and class are the same for 1/0, respectively
		if ( hout[f,"orth"] == 1 && grepl("\\|TE", hout[f,"speciesID"]) ) { 
			if ( gsub(".*\\|", "", hsub[,"V4"]) == gsub(".*\\|", "", msub[,"V6"]) ) { hout[f,"mate_same_te"] <- 1 } else { hout[f,"mate_same_te"] <- 0 } }
	}

	# Now do the same thing but in a mouse-centric manner! Yeeeeeeeeehaw.
	for (f in 1:nrow(mout)) {
		hich <- which(hlo$V7 == mout[f,1]); mich <- which(mlo$V5 == mout[f,1])
		hsub <- hlo[hich,]; msub <- mlo[mich,]

		### ORTHOLOGY
		# If at least 1 msub instance reciprocally maps back to itself with the current mout value, orthologous. Otherwise, not orthologous 
		if ( sum( msub$V5 == msub$V9 ) > 0 ) { mout[f,"orth"] <- 1 
			} else { mout[f,"orth"] <- mout[f,"mate_act"] <- 0; if ( grepl("\\|TE", mout[f,"speciesID"]) ) { mout[f,"mate_same_te"] <- 0 }  }

		### MATE ACTIVITY
		# If not orthologous, 0 (change nothing as done above). If orthologous, check for at least one peak ID in msub or no hich length for 0/1, respectively
		if ( mout[f,"orth"] == 1 && ( sum( msub$V7 != "." ) > 0 || length(hich) > 0 ) ) { mout[f,"mate_act"] <- 1 } else { mout[f,"mate_act"] <- 0 }

		### SAME TE
		# If not orthologous AND TE, 0 (no change). If not TE, -1 (no change). Otherwise, check if the family and class are the same for 1/0, respectively
		if ( mout[f,"orth"] == 1 && grepl("\\|TE", mout[f,"speciesID"]) ) { 
			if ( gsub(".*\\|", "", msub[,"V4"]) == gsub(".*\\|", "", hsub[,"V6"]) ) { mout[f,"mate_same_te"] <- 1 } else { mout[f,"mate_same_te"] <- 0 } }
	}

	# Time to combine the IDs together along with each of the various values to make things easy to see! Also technically not entirely needed 
	# at all, so yay?
	hout$newSpeciesID <- paste0(hout$speciesID, ifelse(hout$orth == 1, "|Orth", "|NonOrth"), ifelse(hout$mate_act == 1, "|Both", "|Self"), ifelse(hout$mate_same_te == 1, "|SameTE", ifelse(hout$mate_same_te == 0, "|DifTE", "|NoTE") ) )
	rownames(hout) <- hout$newSpeciesID

	mout$newSpeciesID <- paste0(mout$speciesID, ifelse(mout$orth == 1, "|Orth", "|NonOrth"), ifelse(mout$mate_act == 1, "|Both", "|Self"), ifelse(mout$mate_same_te == 1, "|SameTE", ifelse(mout$mate_same_te == 0, "|DifTE", "|NoTE") ) )
	rownames(mout) <- mout$newSpeciesID

	hout$sample <- samps[i]; hout$species <- "human"; mout$sample <- samps[i]; mout$species <- "mouse"

	# With everything now done, it's time to combine the various out files to the various win files!
	colnames(hwin) <- c("chr", "slopStart", "slopEnd", "start", "end", "teID", "speciesID", "nearID"); colnames(mwin) <- c("chr", "slopStart", "slopEnd", "start", "end", "teID", "speciesID", "nearID")

	# Instead of merging, add all the different columns to then drop those analyzed from hout/mout into hwin/mwin
	hwin$orth <- mwin$orth <- hwin$mate_act <- mwin$mate_act <- 0; hwin$mate_same_te <- mwin$mate_same_te <- -1 
	hwin$newSpeciesID <- paste0(hwin$speciesID, "|NonOrth|Self", ifelse(grepl("\\|TE", hwin$speciesID), "|difTE", "|NoTE") )
	mwin$newSpeciesID <- paste0(mwin$speciesID, "|NonOrth|Self", ifelse(grepl("\\|TE", mwin$speciesID), "|difTE", "|NoTE") )
	hwin$sample <- samps[i]; hwin$species <- "human"; mwin$sample <- samps[i]; mwin$species <- "mouse"

	# You'll need to iterate through hout/mout to move all the information over. There's probably a better way but eh
	for (f in 1:nrow(hout)) {
		crow <- which(hwin$speciesID == hout[f,"speciesID"]); hwin[crow, colnames(hout)[2:5]] <- hout[f,colnames(hout)[2:5]]
		crow <- which(mwin$speciesID == hout[f,"speciesID"]); mwin[crow, colnames(hout)[2:5]] <- hout[f,colnames(hout)[2:5]]
	}

	for (f in 1:nrow(mout)) {
		crow <- which(hwin$speciesID == mout[f,"speciesID"]); hwin[crow, colnames(mout)[2:5]] <- mout[f,colnames(mout)[2:5]]
		crow <- which(mwin$speciesID == mout[f,"speciesID"]); mwin[crow, colnames(mout)[2:5]] <- mout[f,colnames(mout)[2:5]]
	}

	# Make a lookup table for the speciesID and the newSpeciesID for fast changing of the nearID
	lot <- rbind(hwin[,c("speciesID", "newSpeciesID")], mwin[,c("speciesID", "newSpeciesID")]) %>% distinct(); rownames(lot) <- lot$speciesID

	# Now generate the newNearID based on the newSpeciesID
	hwin$newNearID <- ifelse(hwin$nearID == ".", ".", lot[hwin$nearID,"newSpeciesID"]); mwin$newNearID <- ifelse(mwin$nearID == ".", ".", lot[mwin$nearID,"newSpeciesID"])

	write.table(hwin[,c("chr", "start", "end", "orth", "mate_act", "mate_same_te", "teID", "newSpeciesID", "newNearID")], file = paste0("2_humanTFs/", samps[i], "_5kb_merge"), quote = F, row.names = F, sep = "\t")
	write.table(mwin[,c("chr", "start", "end", "orth", "mate_act", "mate_same_te", "teID", "newSpeciesID", "newNearID")], file = paste0("1_mouseTFs/", samps[i], "_5kb_merge"), quote = F, row.names = F, sep = "\t")

	# Add the syntenic coordinates for each out file
	hint <- read.table(paste0("1_mouseTFs/hg38.mm10.", samps[i], ".int")); mint <- read.table(paste0("2_humanTFs/mm10.hg38.", samps[i], ".int"))
	hint <- hint[order(hint$V5),]; mint <- mint[order(mint$V5),]

	# Remove duplicate speciesID entries and just choose the first one (very infrequent)
	hint <- hint[match(unique(hint$V5), hint$V5),]; mint <- mint[match(unique(mint$V5), mint$V5),]

	rownames(hint) <- hint$V5; rownames(mint) <- mint$V5
	hint$coords <- ifelse(hint$V7 == ".", paste0(hint$V1, ":", hint$V2, "-", hint$V3, "|noPeak"), paste0(hint$V1, ":", hint$V2, "-", hint$V3, "|yesPeak"))
	mint$coords <- ifelse(mint$V7 == ".", paste0(mint$V1, ":", mint$V2, "-", mint$V3, "|noPeak"), paste0(mint$V1, ":", mint$V2, "-", mint$V3, "|yesPeak"))
	hout$syntenicID <- hint[hout$speciesID,"coords"]; mout$syntenicID <- mint[mout$speciesID,"coords"]
	
	rbind(hout, mout)

}, mc.cores = 30)


# Now that the above has been made, combine the rows together and write out the two different files for each species
output <- bind_rows(results)
write.table(subset(output, species == "human")[,c("newSpeciesID", "syntenicID", "sample")], file = paste0(format(Sys.Date(),format="%d%m%y"), "_syntenicPeaks_human.txt"), quote = F, row.names = F, sep = "\t")
write.table(subset(output, species == "mouse")[,c("newSpeciesID", "syntenicID", "sample")], file = paste0(format(Sys.Date(),format="%d%m%y"), "_syntenicPeaks_mouse.txt"), quote = F, row.names = F, sep = "\t")

