#Commands used to compare human and mouse cCREs

#Get human cCRE-TE overlaps and remove any duplicates caused by two (or more) TEs in the same genomic locus
tail -n+2 ../TE_derived_cCREs_human/allHg38ccre_allTE_overlap > allHg38ccre_allTE_overlap
python3 remove_duplicate_overlaps.py allHg38ccre_allTE_overlap allHg38ccre_allTE_overlap_noDups

#Perform liftOver to find syntenic mouse regions to human cCREs (if possible)
cut -f 1-3 allHg38ccre_allTE_overlap_noDups > temp
liftOver -minMatch=0.1 temp hg38ToMm10.over.chain.gz liftOver_hg38_ccre_TE_overlap failed_liftOvers_hg38
#Match human cCREs to mouse syntenic regions
sort -k 1,1 -k2,2n liftOver_hg38_ccre_TE_overlap > liftOver_hg38_ccre_TE_overlap_sorted
bedtools intersect -wao -sorted -a liftOver_hg38_ccre_TE_overlap_sorted -b allMm10ccre_allTE_overlap_noDups > hg38_ccreTE_overlap_mm10_ccreTE
python3 match_elements_after_liftOver_intersect_hg38mm10.py allHg38ccre_allTE_overlap_noDups liftOver_hg38_ccre_TE_overlap failed_liftOvers_hg38 hg38_ccreTE_overlap_mm10_ccreTE liftOver_hg38_ccre_TE_overlap_matched
#Intersect mouse syntenic regions to mouse cCREs and TEs
bedtools intersect -wao -f 0.5 -sorted -a liftOver_hg38_ccre_TE_overlap_sorted -b ccres_mm10_sorted.bed > hg38_ccreTE_overlap_mm10_ccre
bedtools intersect -wao -f 0.5 -sorted -a liftOver_hg38_ccre_TE_overlap_sorted -b <(cut -f 1-6 reAnno_mm10_ENCOTE.bed) > hg38_ccreTE_overlap_mm10_TE
python3 match_elements_after_liftOver_intersect_hg38mm10.py allHg38ccre_allTE_overlap_noDups liftOver_hg38_ccre_TE_overlap failed_liftOvers_hg38 hg38_ccreTE_overlap_mm10_ccre liftOver_hg38_ccre_TE_overlap_matched_mm10ccre
python3 match_elements_after_liftOver_intersect_hg38mm10.py allHg38ccre_allTE_overlap_noDups liftOver_hg38_ccre_TE_overlap failed_liftOvers_hg38 hg38_ccreTE_overlap_mm10_TE liftOver_hg38_ccre_TE_overlap_matched_mm10TE
#Combine info together
python3 parse_liftOver_match_ccre.py liftOver_hg38_ccre_TE_overlap_matched_mm10ccre summary_hg38_mm10_ccre_liftOver
python3 parse_liftOver_match_TE.py liftOver_hg38_ccre_TE_overlap_matched_mm10TE summary_hg38_mm10_TE_liftOver
python3 quantify_human_mouse_ccre_TEs_v2.py summary_hg38_mm10_ccre_liftOver summary_hg38_mm10_TE_liftOver ../age_calculation/kimura_distance_hg38.tsv ../age_calculation/kimura_distance_mm10.tsv reAnno_hg38.bed final_info_human > final_summary_human 2> final_info_human_nonMatching_TE_cCREs
python3 edit_summary_forR.py final_summary_human forR_final_summary_human
cat <(tail -n18 final_summary_human | head -n6 | sed 's/#//g' | cut -f 1-4,6) <(tail -n10 final_summary_human | head -n5 | sed 's/#//g') > forR_summary_novel_human
#Add clade info
python3 assign_TE_clade.py final_info_human /scratch/adu/encode_te_analysis/dfam/shared_human_rodent_TEs_hg38 /scratch/adu/encode_te_analysis/dfam/primate_human_lineage_TEs final_info_human_wClade

#Do reciprocal analysis starting from mouse
mkdir reciprocal_mouse
bedtools intersect -wao -f 0.5 -sorted -a ccres_mm10_sorted.bed -b <(cut -f 1-7 reAnno_mm10_ENCOTE.bed) > reciprocal_mouse/allMm10ccre_allTE_overlap
python3 remove_duplicate_overlaps.py reciprocal_mouse/allMm10ccre_allTE_overlap reciprocal_mouse/allMm10ccre_allTE_overlap_noDups

cut -f 1-3 reciprocal_mouse/allMm10ccre_allTE_overlap_noDups > temp
liftOver -minMatch=0.1 temp mm10ToHg38.over.chain.gz reciprocal_mouse/liftOver_mm10_ccre_TE_overlap reciprocal_mouse/failed_liftOvers_mm10

sort -k 1,1 -k2,2n reciprocal_mouse/liftOver_mm10_ccre_TE_overlap > reciprocal_mouse/liftOver_mm10_ccre_TE_overlap_sorted
bedtools intersect -wao -f 0.5 -sorted -a reciprocal_mouse/liftOver_mm10_ccre_TE_overlap_sorted -b ccres_hg38_sorted.bed > reciprocal_mouse/mm10_ccreTE_overlap_hg38_ccre
bedtools intersect -wao -f 0.5 -sorted -a reciprocal_mouse/liftOver_mm10_ccre_TE_overlap_sorted -b <(cut -f 1-6 reAnno_hg38.bed) > reciprocal_mouse/mm10_ccreTE_overlap_hg38_TE
python3 match_elements_after_liftOver_intersect_hg38mm10.py reciprocal_mouse/allMm10ccre_allTE_overlap_noDups reciprocal_mouse/liftOver_mm10_ccre_TE_overlap reciprocal_mouse/failed_liftOvers_mm10 reciprocal_mouse/mm10_ccreTE_overlap_hg38_ccre reciprocal_mouse/liftOver_mm10_ccre_TE_overlap_matched_hg38_ccre
python3 match_elements_after_liftOver_intersect_hg38mm10.py reciprocal_mouse/allMm10ccre_allTE_overlap_noDups reciprocal_mouse/liftOver_mm10_ccre_TE_overlap reciprocal_mouse/failed_liftOvers_mm10 reciprocal_mouse/mm10_ccreTE_overlap_hg38_TE reciprocal_mouse/liftOver_mm10_ccre_TE_overlap_matched_hg38_TE

python3 parse_liftOver_match_ccre.py reciprocal_mouse/liftOver_mm10_ccre_TE_overlap_matched_hg38_ccre reciprocal_mouse/summary_mm10_hg38_ccre_liftOver
python3 parse_liftOver_match_TE.py reciprocal_mouse/liftOver_mm10_ccre_TE_overlap_matched_hg38_TE reciprocal_mouse/summary_mm10_hg38_TE_liftOver
python3 quantify_mouse_human_ccre_TEs_v2.py reciprocal_mouse/summary_mm10_hg38_ccre_liftOver reciprocal_mouse/summary_mm10_hg38_TE_liftOver ../age_calculation/kimura_distance_hg38.tsv ../age_calculation/kimura_distance_mm10.tsv reAnno_mm10_ENCOTE.bed reciprocal_mouse/final_info_mouse > reciprocal_mouse/final_summary_mouse 2> reciprocal_mouse/final_info_mouse_nonMatching_TE_cCREs
python3 edit_summary_forR.py reciprocal_mouse/final_summary_mouse reciprocal_mouse/forR_final_summary_mouse
cat <(tail -n18 reciprocal_mouse/final_summary_mouse | head -n6 | sed 's/#//g' | cut -f 1-4,6) <(tail -n10 reciprocal_mouse/final_summary_mouse | head -n5 | sed 's/#//g') > reciprocal_mouse/forR_summary_novel_mouse

#phastCons100way (and phyloP100way) analysis
mkdir phastCons_phylop
cd phastCons_phylop
ln -s ../*.py .
#Start by making bedGraph from bigWig file
bigWigToBedGraph hg38.phastCons100way.bw hg38.phastCons100way.bedGraph
bigWigToBedGraph hg38.phyloP100way.bw hg38.phyloP100way.bedGraph
mkdir sorting_files
mkdir sorting_files/temp
grep "chr1" hg38.phastCons100way.bedGraph > sorting_files/chr1s_phastCons_sorted.bedGraph
grep "chr2" hg38.phastCons100way.bedGraph > sorting_files/chr2s_phastCons_sorted.bedGraph
grep -v "chr1" hg38.phastCons100way.bedGraph | grep -v "chr2" > sorting_files/restChrs_phastCons_sorted.bedGraph
sort --parallel=12 -k1,1 -k2,2n -T sorting_files/temp/ sorting_files/chr1s_phastCons_sorted.bedGraph > hg38.phastCons100way.sorted.bedGraph
sort --parallel=12 -k1,1 -k2,2n -T sorting_files/temp/ sorting_files/chr2s_phastCons_sorted.bedGraph >> hg38.phastCons100way.sorted.bedGraph
sort --parallel=12 -k1,1 -k2,2n -T sorting_files/temp/ sorting_files/restChrs_phastCons_sorted.bedGraph >> hg38.phastCons100way.sorted.bedGraph
grep "chr1" hg38.phyloP100way.bedGraph > sorting_files/chr1s_phyloP_sorted.bedGraph
grep "chr2" hg38.phyloP100way.bedGraph > sorting_files/chr2s_phyloP_sorted.bedGraph
grep -v "chr1" hg38.phyloP100way.bedGraph | grep -v "chr2" > sorting_files/restChrs_phyloP_sorted.bedGraph
sort --parallel=12 -k1,1 -k2,2n -T sorting_files/temp/ sorting_files/chr1s_phyloP_sorted.bedGraph > hg38.phyloP100way.sorted.bedGraph
sort --parallel=12 -k1,1 -k2,2n -T sorting_files/temp/ sorting_files/chr2s_phyloP_sorted.bedGraph >> hg38.phyloP100way.sorted.bedGraph
sort --parallel=12 -k1,1 -k2,2n -T sorting_files/temp/ sorting_files/restChrs_phyloP_sorted.bedGraph >> hg38.phyloP100way.sorted.bedGraph
rm hg38.phastCons100way.bedGraph
rm hg38.phyloP100way.bedGraph
rm -r sorting_files
#Get mean phastCons and phyloP 100way scores for cCREs
bedtools intersect -sorted -wao -a hg38_ccres_sorted.bed -b hg38.phastCons100way.sorted.bedGraph > phastCons100way_hg38_ccres
python3 combine_phylop_phastcons_scores_ccres.py phastCons100way_hg38_ccres phastCons100way_hg38_ccres_combined > noScores_phastCons100way_hg38_ccres
bedtools intersect -sorted -wao -a hg38_ccres_sorted.bed -b hg38.phyloP100way.sorted.bedGraph > phyloP100way_hg38_ccres
python3 combine_phylop_phastcons_scores_ccres.py phyloP100way_hg38_ccres phyloP100way_hg38_ccres_combined > noScores_phyloP100way_hg38_ccres
#Get mean phastCons and phyloP 100way scores for TEs
bedtools intersect -sorted -wao -a reAnno_hg38.bed -b hg38.phastCons100way.sorted.bedGraph > phastCons100way_hg38_TEs
python3 combine_phylop_phastcons_scores_TEs.py phastCons100way_hg38_TEs phastCons100way_hg38_TEs_combined > noScores_phastCons100way_hg38_TEs
bedtools intersect -sorted -wao -a reAnno_hg38.bed -b hg38.phyloP100way.sorted.bedGraph > phyloP100way_hg38_TEs
python3 combine_phylop_phastcons_scores_TEs.py phyloP100way_hg38_TEs phyloP100way_hg38_TEs_combined > noScores_phyloP100way_hg38_TEs

#Assign phastCons and phyloP 100way scores to cCREs
python3 assign_phylop_phastcons_scores_ccres.py ../final_info_human_wClade phyloP100way_hg38_ccres_combined phastCons100way_hg38_ccres_combined final_info_human_wClade_100way
python3 assign_phylop_phastcons_scores_ccres.py ../final_info_human phyloP100way_hg38_ccres_combined phastCons100way_hg38_ccres_combined final_info_human_100way
#Remove CTCF-bound annotation for R
sed 's/,CTCF-bound//g' final_info_human_wClade_100way > forR_final_info_human_wClade_100way

#Compare phastCons scores between TEs based on synteny and conserved cCRE annotation
python3 calc_phylop_phastcons_perSubfamily_compare_orthologTEs.py ../final_info_human_wClade phyloP100way_hg38_TEs_combined phastCons100way_hg38_TEs_combined compare_orthologousTE_phyloP_phastCons_100way > compare_orthologousTE_phyloP_phastCons_100way_allValues
python3 calc_phylop_phastcons_perSubfamily_compare_orthologTEs_cCRE_vs_non.py ../final_info_human_wClade phyloP100way_hg38_TEs_combined phastCons100way_hg38_TEs_combined compare_orthologousTE_vs_non_phyloP_phastCons_100way > compare_orthologousTE_vs_non_phyloP_phastCons_100way_allValues

