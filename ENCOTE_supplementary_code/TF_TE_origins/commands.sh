#Commands used to estimate TF motif origins

mkdir forR
ln -s ../TE_derived_cCREs_human/separate_TEsubfamily_bed/ .
ln -s ../TE_derived_cCREs_human/TEsubfamilies_grouped/ .
ln -s ../reference_files/list_hg38_TE_subfamilies .
ln -s ../reference_files/list_TE_classes_human .
ln -s ../reference_files/archetype_clusters_info.txt .
ln -s ../reference_files/archetype_motif_info.txt .
ln -s ../reference_files/kimura_distance_hg38.tsv .
ln -s ../reference_files/hg38_ccres_sorted.bed .
ln -s ../reference_files/final_hg38_consensus_seq.fa hg38TE.consensus.fa
ln -s ../reference_files/alignFastas/ .

grep ">" hg38TE.consensus.fa | sed 's/>//g' | sort > hg38TE.consensus.txt

#Scan for TF motifs in TE consensus sequences
fimo --oc motif_scan_consensus --verbosity 1 --thresh 1E-4 /opt/apps/meme/5.4.1/db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme hg38TE.consensus.fa
python3 parse_consensus_motif_scan.py motif_scan_consensus/fimo.tsv list_hg38_TE_subfamilies consensus_motifs.txt

#Get sequences of TEs separated by hg38 cCRE overlap annotation for each subfamily
mkdir overlaps_TEsubfamilies_hg38_ccres
mkdir sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap
for i in `cat list_hg38_TE_subfamilies`; do
        bedtools intersect -sorted -u -F 0.5 -a separate_TEsubfamily_bed/"$i".bed -b hg38_ccres_sorted.bed > overlaps_TEsubfamilies_hg38_ccres/"$i".overlap &
        bedtools intersect -sorted -v -F 0.5 -a separate_TEsubfamily_bed/"$i".bed -b hg38_ccres_sorted.bed > overlaps_TEsubfamilies_hg38_ccres/"$i".nonOverlap &
        wait
done
for group in `ls TEsubfamilies_grouped/`; do
        for i in `cat TEsubfamilies_grouped/$group`; do
                fastaFromBed -s -bed overlaps_TEsubfamilies_hg38_ccres/"$i".overlap -fi ../reference_files/hg38.fa -fo sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_overlap.fa &
                fastaFromBed -s -bed overlaps_TEsubfamilies_hg38_ccres/"$i".nonOverlap -fi ../reference_files/hg38.fa -fo sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_nonOverlap.fa &
        done
        wait
done

#First run AME motif enrichment for cCRE overlapping (foreground) TEs vs. (background) non-overlapping TEs
mkdir AME_enriched_motifs_hg38_ccres
for group in `ls TEsubfamilies_grouped/`; do
        for i in `cat TEsubfamilies_grouped/$group`; do
                ame --oc AME_enriched_motifs_hg38_ccres/"$i" --control sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_nonOverlap.fa sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_overlap.fa /meme/5.4.1/db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme &
        done
        wait
done
#Get more stringent criteria for motif enrichment
#At least 50% of motifs from a motif archetype/cluster are enriched in a subfamily
mkdir stringent_enriched_motifs
for i in `cat list_TE_subfamilies`; do
        python3 get_stringent_enriched_motifs.py AME_enriched_motifs_hg38_ccres/"$i"/ame.tsv archetype_motif_info.txt archetype_clusters_info.txt stringent_enriched_motifs/"$i".enriched_motifs
done
#Scan all TEs of a subfamily for enriched motifs
mkdir motif_scan_individual_TEs
mkdir sequences_TEsubfamilies
#counter used to wait on every 12th fimo command
counter=0
for group in `ls TEsubfamilies_grouped/`; do
        for i in `cat TEsubfamilies_grouped/$group`; do
                mkdir motif_scan_individual_TEs/"$i"
                cat sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_overlap.fa sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_nonOverlap.fa > sequences_TEsubfamilies/"$i".fa
                for motif in `cut -f 3 stringent_enriched_motifs/"$i".enriched_motifs | tail -n+2`; do
                        ((counter+=1))
                        fimo --oc motif_scan_individual_TEs/"$i"/"$motif" --max-stored-scores 1000000 --verbosity 1 --thresh 1E-4 --motif "$motif" /meme/5.4.1/db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme sequences_TEsubfamilies/"$i".fa &
                        if [ $counter == 12 ]; then
                                wait
                                counter=0
                        fi
                done
        done
done
#Further filter for enriched motifs that are significant in fisher's exact test using FIMO identified motifs
mkdir annotations_enriched_motifs_perSubfamily
for group in `ls TEsubfamilies_grouped/`; do
        for i in `cat TEsubfamilies_grouped/$group`; do
                python3 get_enriched_motif_annotations.py separate_TEsubfamily_bed/"$i".bed motif_scan_individual_TEs/"$i"/ annotations_enriched_motifs_perSubfamily/"$i".motif_annotations &
        done
        wait
done
mkdir fishers_cCRE_enriched_motifs_perSubfamily
for group in `ls TEsubfamilies_grouped/`; do
        for i in `cat TEsubfamilies_grouped/$group`; do
                python3 find_motifs_assoc_wIntersect_fishers.py annotations_enriched_motifs_perSubfamily/"$i".motif_annotations overlaps_TEsubfamilies_hg38_ccres/"$i".overlap overlaps_TEsubfamilies_hg38_ccres/"$i".nonOverlap fishers_cCRE_enriched_motifs_perSubfamily/"$i".fishers &
        done
        wait
done
#Get high confidence enriched motifs from Fisher's exact test results
mkdir high_confidence_enriched_motifs
for i in `cat list_TE_subfamilies`; do
        python3 get_high_confident_enriched_motifs.py stringent_enriched_motifs/"$i".enriched_motifs fishers_cCRE_enriched_motifs_perSubfamily/"$i".fishers high_confidence_enriched_motifs/"$i".confident_motifs
done

#Next, find TE subfamilies that have significant length distribution difference between foreground and background elements
mkdir KStests_length_distribution
for i in `cat list_TE_subfamilies`; do
        python3 get_significant_coverage_subfamilies_kstest.py sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_overlap.fa sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_nonOverlap.fa "$i" >> KStests_length_distribution/all_background_elements 2>> KStests_length_distribution/errors_all_background.txt
done
python3 get_length_biased_subfamilies_kstest_noCorrect.py KStests_length_distribution/all_background_elements list_length_biased_subfamilies_kstest_conservative > list_length_nonbiased_subfamilies_kstest_conservative
#Get list of TE subfamilies that can get enough background elements after matching length distribution of foreground (don't consider TE subfamilies that can't possibly have matched background)
mkdir random_background_cCRE_elements_max_v3
for i in `cat list_length_biased_subfamilies_kstest_conservative`; do
        mkdir random_background_cCRE_elements_max_v3/"$i"
        python3 randomize_background_cCRE_motif_enrichment_v3.py sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_overlap.fa sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_nonOverlap.fa 10 random_background_cCRE_elements_max_v3/"$i"/"$i".random 2>> random_background_cCRE_elements_max_v3/errors_random_background.txt
done
grep "Foreground" random_background_cCRE_elements_max_v3/errors_random_background.txt | cut -d" " -f 2 | cut -f 1 | cut -d"/" -f 2 | cut -d"_" -f 1 > random_background_cCRE_elements_max_v3/list_filtered_subfamilies
cut -d" " -f 1 random_background_cCRE_elements_max_v3/errors_random_background.txt | cut -f 1 | sort | uniq | grep -v "scipy" | grep -v "Foreground" | grep -v "No" | tail -n+2 >> random_background_cCRE_elements_max_v3/list_filtered_subfamilies
comm -23 <(sort list_length_biased_subfamilies_kstest_conservative ) random_background_cCRE_elements_max_v3/list_filtered_subfamilies > list_length_biased_subfamilies_kstest_conservative_wRandom_v3
#Do motif enrichment analysis for length controlled background elements
mkdir random_background_cCRE_elements_overlaps_max_v3
for j in `seq 1 10`; do
        mkdir random_background_cCRE_elements_overlaps_max_v3/random"$j"
        cd random_background_cCRE_elements_overlaps_max_v3/random"$j"/
        for i in `ls ../../overlaps_TEsubfamilies_hg38_ccres/*.overlap | cut -d"/" -f 4`; do
                ln -s ../../overlaps_TEsubfamilies_hg38_ccres/"$i" .
        done
        cd ../..
        for i in `cat list_length_biased_subfamilies_kstest_conservative_wRandom_v3`; do
                python3 transform_randomBackground_fasta_toOverlap.py random_background_cCRE_elements_max_v3/"$i"/"$i".random"$j" random_background_cCRE_elements_overlaps_max_v3/random"$j"/"$i".nonOverlap
        done
done
#Find enriched motifs using AME using random background
mkdir AME_enriched_motifs_hg38_ccres_randomBackground_max_v3
for i in `cat list_length_biased_subfamilies_kstest_conservative_wRandom_v3`; do
        mkdir AME_enriched_motifs_hg38_ccres_randomBackground_max_v3/"$i"
        for j in `seq 1 10`; do
                ame --oc AME_enriched_motifs_hg38_ccres_randomBackground_max_v3/"$i"/random"$j" --control random_background_cCRE_elements_max_v3/"$i"/"$i".random"$j" sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_overlap.fa /meme/5.4.1/db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme &
        done
        wait
done
wait
#Get stringent set of motifs for AME enriched motifs
mkdir stringent_enriched_motifs_randomBackground_max_v3
for i in `cat list_length_biased_subfamilies_kstest_conservative_wRandom_v3`; do
        python3 get_stringent_enriched_motifs_randomBackground_v2.py AME_enriched_motifs_hg38_ccres_randomBackground_max_v3/"$i"/ archetype_motif_info.txt archetype_clusters_info.txt stringent_enriched_motifs_randomBackground_max_v3/"$i".enriched_motifs
done
#Scan individual elements for enriched motifs
mkdir motif_scan_individual_TEs_randomMax_v3
counter=0
for i in `cat list_length_biased_subfamilies_kstest_conservative_wRandom_v3`; do
        mkdir motif_scan_individual_TEs_randomMax_v3/"$i"
        for motif in `cut -f 3 stringent_enriched_motifs_randomBackground_max_v3/"$i".enriched_motifs | tail -n+2`; do
                OUT_DIR=motif_scan_individual_TEs_randomMax_v3/"$i"/"$motif"
                if test -e "$OUT_DIR"; then
                        continue
                else
                        ((counter+=1))
                        fimo --o motif_scan_individual_TEs_randomMax_v3/"$i"/"$motif" --max-stored-scores 1000000 --verbosity 1 --thresh 1E-4 --motif "$motif" /meme/5.4.1/db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme sequences_TEsubfamilies/"$i".fa &
                        if [ $counter == 12 ]; then
                                wait
                        fi
                fi
        done
done
wait
echo "Finished individual element motif scans"
#Annotate elements for motif presence and perform Fisher's exact test
mkdir annotations_enriched_motifs_perSubfamily_randomBackground_max_v3
for i in `cat list_length_biased_subfamilies_kstest_conservative_wRandom_v3`; do
        python3 get_enriched_motif_annotations.py separate_TEsubfamily_bed/"$i".bed motif_scan_individual_TEs_randomMax_v3/"$i"/ annotations_enriched_motifs_perSubfamily_randomBackground_max_v3/"$i".motif_annotations
done
mkdir fishers_cCRE_enriched_motifs_perSubfamily_randomBackground_max_v3
for i in `cat list_length_biased_subfamilies_kstest_conservative_wRandom_v3`; do
        mkdir fishers_cCRE_enriched_motifs_perSubfamily_randomBackground_max_v3/"$i"
        for j in `seq 1 10`; do
                python3 find_motifs_assoc_wIntersect_fishers_randomBackground.py stringent_enriched_motifs_randomBackground_max_v3/"$i".enriched_motifs annotations_enriched_motifs_perSubfamily_randomBackground_max_v3/"$i".motif_annotations overlaps_TEsubfamilies_hg38_ccres/"$i".overlap random_background_cCRE_elements_overlaps_max_v3/random"$j"/"$i".nonOverlap fishers_cCRE_enriched_motifs_perSubfamily_randomBackground_max_v3/"$i"/"$i".random"$j".fishers >> fishers_cCRE_enriched_motifs_perSubfamily_randomBackground_max_v3/errors_fishers_randomBackground.txt &
        done
        wait
done
wait
#Get high confident motifs from random background
mkdir high_confidence_enriched_motifs_randomBackground_max_v3
for i in `cat list_length_biased_subfamilies_kstest_conservative_wRandom_v3`; do
        python3 get_high_confident_enriched_motifs_randomBackground.py stringent_enriched_motifs_randomBackground_max_v3/"$i".enriched_motifs fishers_cCRE_enriched_motifs_perSubfamily_randomBackground_max_v3/"$i"/ high_confidence_enriched_motifs_randomBackground_max_v3/"$i".confident_motifs
done
#Combine high confident motifs from random background and full background subfamilies
mkdir combined_high_confidence_enriched_motifs_v3
cd combined_high_confidence_enriched_motifs_v3
for i in `cat ../list_length_nonbiased_subfamilies_kstest_conservative`; do
        ln -s ../high_confidence_enriched_motifs/"$i".confident_motifs .
done
for i in `cat ../list_length_biased_subfamilies_kstest_conservative_wRandom_v3`; do
        ln -s ../high_confidence_enriched_motifs_randomBackground_max_v3/"$i".confident_motifs .
done
cd ..

#Next, control for consensus coverage between foreground and background (to increase specificity of identified cCRE associated TF motifs)
#Find TE subfamilies that can have similar foreground and background consensus coverage
mkdir controlled_coverage_background
mkdir controlled_coverage_background_forKS
for i in `cat list_hg38_TE_subfamilies`; do
        OUT_FILE=controlled_coverage_background/"$i"_background.fa
        if test -f "$OUT_FILE"; then
                continue
        else
                python3 select_coverage_background_cCRE_motif_enrichment_v4.py alignFastas/"$i".alignFasta overlaps_TEsubfamilies_hg38_ccres/"$i".overlap overlaps_TEsubfamilies_hg38_ccres/"$i".nonOverlap controlled_coverage_background/"$i" > controlled_coverage_background_forKS/"$i" 2>> log_controlled_coverage_background
        fi
done
comm -23 <(sort list_hg38_TE_subfamilies) <(cut -f 1 log_controlled_coverage_background | sort) > list_coverage_controlled_subfamilies
#Find enriched motifs using AME using coverage controlled background elements
mkdir AME_enriched_motifs_hg38_ccres_coverageControl
counter=0
for i in `cat list_coverage_controlled_subfamilies`; do
        ((counter+=1))
        ame --oc AME_enriched_motifs_hg38_ccres_coverageControl/"$i" --control controlled_coverage_background/"$i"_background.fa sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_overlap.fa /meme/5.4.1/db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme &
        if [ $counter == 12 ]; then
                wait
        fi
done
wait
#Get stringent set of motifs for AME enriched motifs
mkdir stringent_enriched_motifs_coverageControl
for i in `cat list_coverage_controlled_subfamilies`; do
        python3 get_stringent_enriched_motifs.py AME_enriched_motifs_hg38_ccres_coverageControl/"$i"/ame.tsv archetype_motif_info.txt archetype_clusters_info.txt stringent_enriched_motifs_coverageControl/"$i".enriched_motifs
done
#Scan individual elements for enriched motifs (that were not previously already done)
counter=0
for i in `cat list_coverage_controlled_subfamilies`; do
        OUT_DIR=motif_scan_individual_TEs_randomMax_v3/"$i"/
        if ! test -e "$OUT_DIR"; then
                mkdir motif_scan_individual_TEs_randomMax_v3/"$i"
        fi
        for motif in `cut -f 3 stringent_enriched_motifs_coverageControl/"$i".enriched_motifs | tail -n+2`; do
                OUT_DIR=motif_scan_individual_TEs_randomMax_v3/"$i"/"$motif"
                if test -e "$OUT_DIR"; then
                        continue
                else
                        ((counter+=1))
                        fimo --o motif_scan_individual_TEs_randomMax_v3/"$i"/"$motif" --max-stored-scores 1000000 --verbosity 1 --thresh 1E-4 --motif "$motif" /meme/5.4.1/db/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme sequences_TEsubfamilies/"$i".fa &
                        if [ $counter == 12 ]; then
                                wait
                        fi
                fi
        done
done
wait
echo "Finished individual element motif scans"
#Perform Fisher's exact test using previously annotated motifs
mkdir fishers_cCRE_enriched_motifs_perSubfamily_coverageControl
for i in `cat list_coverage_controlled_subfamilies`; do
        IN_FILE1=stringent_enriched_motifs_randomBackground_max_v3/"$i".enriched_motifs
        IN_FILE2=stringent_enriched_motifs/"$i".enriched_motifs
        if test -e "$IN_FILE1"; then
                python3 find_motifs_assoc_wIntersect_fishers_v3.py stringent_enriched_motifs_randomBackground_max_v3/"$i".enriched_motifs annotations_enriched_motifs_perSubfamily_randomBackground_max_v3/"$i".motif_annotations overlaps_TEsubfamilies_hg38_ccres/"$i".overlap controlled_coverage_background/"$i"_background.bed fishers_cCRE_enriched_motifs_perSubfamily_coverageControl/"$i".fishers >> fishers_cCRE_enriched_motifs_perSubfamily_coverageControl/errors_fishers.txt
        elif test -e "$IN_FILE2"; then
                python3 find_motifs_assoc_wIntersect_fishers_v3.py stringent_enriched_motifs/"$i".enriched_motifs annotations_enriched_motifs_perSubfamily/"$i".motif_annotations overlaps_TEsubfamilies_hg38_ccres/"$i".overlap controlled_coverage_background/"$i"_background.bed fishers_cCRE_enriched_motifs_perSubfamily_coverageControl/"$i".fishers >> fishers_cCRE_enriched_motifs_perSubfamily_coverageControl/errors_fishers.txt
        else
                continue
        fi
done
#Get high confident motifs that length controlled random background and consensus coverage controlled background agree on
mkdir combined_high_confidence_enriched_motifs_coverageControl
for i in `cat list_coverage_controlled_subfamilies`; do
        IN_FILE1=fishers_cCRE_enriched_motifs_perSubfamily_coverageControl/"$i".fishers
        IN_FILE2=combined_high_confidence_enriched_motifs_v3/"$i".confident_motifs
        if ! test -e "$IN_FILE1"; then
                continue
        elif ! test -e "$IN_FILE2"; then
                continue
        fi
        python3 get_shared_high_confident_motifs.py combined_high_confidence_enriched_motifs_v3/"$i".confident_motifs fishers_cCRE_enriched_motifs_perSubfamily_coverageControl/"$i".fishers combined_high_confidence_enriched_motifs_coverageControl/"$i".confident_motifs
done

#Lastly, calculate ancestral origin percentage for all high confident motifs from length and coverage controlled background
mkdir final_ancOrigin_calculation_coverageControl
python3 calc_confident_motif_origin_consensus_regionOnly_v3.py combined_high_confidence_enriched_motifs_coverageControl/ consensus_motifs.txt hg38TE.consensus.txt list_TE_classes_human kimura_distance_hg38.tsv individual_TE_consensus_motif_location_info/ final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perSubfamily > final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perMotif 2> final_ancOrigin_calculation_coverageControl/errors_regionOnly_confident.txt
sed -e 's/\-int\tLTR/\-int\tERV-int/g' -e 's/MamGypsy2\-I\tLTR/MamGypsy2\-I\tERV-int/g' final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perSubfamily > final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perSubfamily_sepInts
python3 match_subfamily_info.py final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perSubfamily_sepInts length_human_TEsubfamily_consensus final_ancOrigin_calculation_coverageControl/regionOnly_length_human_TEsubfamily_confidentMotifs_consensus_originAnc
#Compare to random motif selection
python3 randomize_enrichedMotif_ancOrigin_test_regionOnly.py final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perSubfamily_sepInts consensus_motifs.txt archetype_motif_info.txt individual_TE_consensus_motif_location_info/ final_ancOrigin_calculation_coverageControl/regionOnly_random_confidentMotifs_consensus_perSubfamily > final_ancOrigin_calculation_coverageControl/regionOnly_random_confidentMotifs_consensus_perRandom
paste final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perSubfamily_sepInts <(cut -f 6,7 final_ancOrigin_calculation_coverageControl/regionOnly_random_confidentMotifs_consensus_perSubfamily ) > final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perSubfamily_sepInts_compareRandom
#Get the number/percentage of motif archetypes that are not present in the consensus sequence (i.e. they are more likely to be completely novel)
python3 get_nonconsensus_enriched_motifs.py final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perMotif kimura_distance_hg38.tsv final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perMotif_nonConsensus
sed -i -e 's/\-int\tLTR/\-int\tERV-int/g' -e 's/MamGypsy2\-I\tLTR/MamGypsy2\-I\tERV-int/g' final_ancOrigin_calculation_coverageControl/regionOnly_originAnc_confidentMotifs_consensus_perMotif_nonConsensus

#Get percentage of elements with cCRE enriched motifs
python3 get_percent_elements_wMotifs_coverageControl.py combined_high_confidence_enriched_motifs_coverageControl kimura_distance_hg38.tsv length_human_TEsubfamily_consensus final_ancOrigin_calculation_coverageControl/percent_motifs
python3 get_mean_percent_elements_wMotifs.py final_ancOrigin_calculation_coverageControl/percent_motifs_coverageControl final_ancOrigin_calculation_coverageControl/percent_motifs_coverageControl_perSubfamily
python3 get_mean_percent_elements_wMotifs.py final_ancOrigin_calculation_coverageControl/percent_motifs_lengthControl final_ancOrigin_calculation_coverageControl/percent_motifs_lengthControl_perSubfamily

#Get consensus coverage for TE subfamilies when not directly controlling for consensus coverage
#First get consensus coverage for all TEs
mkdir coverage_foreground_cCREs
python3 calc_foreground_distribution.py hg38TE.consensus.txt consensus_motifs.txt overlaps_TEsubfamilies_hg38_ccres/ alignFastas/ coverage_foreground_cCREs/cCRE_coverage > exclude_motifs_foreground_coverage
#Next get consensus coverage for random background TEs that match length distribution of foreground
#Randomly select background TEs so that the length distribution matches foreground TEs
mkdir random_background_cCRE_elements_max_v3
for i in `cat list_length_biased_subfamilies_kstest_conservative`; do
        mkdir random_background_cCRE_elements_max_v3/"$i"
        python3 randomize_background_cCRE_motif_enrichment_v3.py sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_overlap.fa sequences_TEsubfamilies_separatedBy_hg38_ccre_overlap/"$i"_nonOverlap.fa 10 random_background_cCRE_elements_max_v3/"$i"/"$i".random 2>> random_background_cCRE_elements_max_v3/errors_random_background.txt
done
grep "Foreground" random_background_cCRE_elements_max_v3/errors_random_background.txt | cut -d" " -f 2 | cut -f 1 | cut -d"/" -f 2 | cut -d"_" -f 1 > random_background_cCRE_elements_max_v3/list_filtered_subfamilies
cut -d" " -f 1 random_background_cCRE_elements_max_v3/errors_random_background.txt | cut -f 1 | sort | uniq | grep -v "scipy" | grep -v "Foreground" | grep -v "No" | tail -n+2 >> random_background_cCRE_elements_max_v3/list_filtered_subfamilies
comm -23 <(sort list_length_biased_subfamilies_kstest_conservative ) random_background_cCRE_elements_max_v3/list_filtered_subfamilies > list_length_biased_subfamilies_kstest_conservative_wRandom_v3
#Calculate the number of TE copies with bases at each position for foreground and background, then calculate relative enrichment at each position for foreground over background
mkdir random_background_cCRE_elements_coverage_max_v3
for j in `seq 1 10`; do
        mkdir random_background_cCRE_elements_coverage_max_v3/random"$j"
        python3 calc_foreground_distribution_wLenThresh_randomBackground.py list_length_biased_subfamilies_kstest_conservative_wRandom_v3 consensus_motifs.txt random_background_cCRE_elements_overlaps_max_v3/random"$j"/ alignFastas/ random_background_cCRE_elements_coverage_max_v3/random"$j"/cCRE_coverage > random_background_cCRE_elements_coverage_max_v3/exclude_motifs_foreground_coverage"$j"
done
#Get consensus coverage for TEs
mkdir length_controlled_random_background_coverage_max_v3
cd length_controlled_random_background_coverage_max_v3/
for i in `seq 1 10`; do
        mkdir random"$i"
        cd random"$i"/
        for subfamily in `cat ../list_hg38_TE_subfamilies`; do
                if test -e ../random_background_cCRE_elements_coverage_max_v3/random"$i"/cCRE_coverage_"$subfamily"; then
                        ln -s ../random_background_cCRE_elements_coverage_max_v3/random"$i"/cCRE_coverage_"$subfamily" .
                elif test -e ../coverage_foreground_cCREs/cCRE_coverage_"$subfamily"; then
                        ln -s ../coverage_foreground_cCREs/cCRE_coverage_"$subfamily" .
                else
                        continue
                fi
        done
        cd ..
done
cd ..
#Get consensus coverage enrichment averages (mean) in different TE classes for all TE subfamilies
for i in `seq 1 10`; do
        python3 summarize_consensus_coverage_inputSubfamilies_byClass.py length_controlled_random_background_coverage_max_v3/random"$i"/ list_TE_subfamilies list_TE_classes length_controlled_random_background_coverage_max_v3/TEclass_consensus_coverage_normEnrichment_"$i"_Class 2> length_controlled_random_background_coverage_max_v3/problems_consensus_coverage_byClass"$i"
        python3 summarize_consensus_coverage_inputSubfamilies_byFamily.py length_controlled_random_background_coverage_max_v3/random"$i"/ list_TE_subfamilies list_TE_classes length_controlled_random_background_coverage_max_v3/TEfamily_consensus_coverage_normEnrichment"$i" 2> length_controlled_random_background_coverage_max_v3/problems_consensus_coverage_byFamily"$i"
done
python3 consolidate_normEnrichment_consensus_coverage.py length_controlled_random_background_coverage_max_v3/TEclass_consensus_coverage_normEnrichment*_Class length_controlled_random_background_coverage_max_v3/combined_normEnrichment_TEclass
python3 consolidate_normEnrichment_consensus_coverage.py length_controlled_random_background_coverage_max_v3/TEfamily_consensus_coverage_normEnrichment*_FamilyOnly length_controlled_random_background_coverage_max_v3/combined_normEnrichment_TEfamily
python3 consolidate_normEnrichment_consensus_coverage.py length_controlled_random_background_coverage_max_v3/TEfamily_consensus_coverage_normEnrichment*_FamilyClass length_controlled_random_background_coverage_max_v3/combined_normEnrichment_TEfamily_wClass

#Find where cCREs overlap in TE families
counter=0
for i in `cat list_hg38_TE_subfamilies`; do
        ((counter+=1))
        bedtools intersect -sorted -wo -F 0.5 -a separate_TEsubfamily_bed/"$i".bed -b hg38_ccres_sorted.bed > overlaps_TEsubfamilies_hg38_ccres/"$i".region_wInfo &
        if [ $counter == 10 ]; then
                counter=0
                wait
        fi
done
mkdir ccreOnly_overlap_coverage
python3 calc_foreground_distribution_regionOnly_overlap.py list_TE_subfamilies overlaps_TEsubfamilies_hg38_ccres/ alignFastas/ ccreOnly_overlap_coverage/coverage > ccreOnly_overlap_coverage/KStests_ccreOnly_overlap_vs_full 2> ccreOnly_overlap_coverage/errors_coverage
python3 consolidate_ccreOnly_overlap_coverage.py ccreOnly_overlap_coverage/ list_hg38_TE_subfamilies list_TE_classes_human ccreOnly_overlap_coverage/TEfamily_ccreOnly_overlap_coverage

