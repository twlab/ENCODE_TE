#Commands used to quantify human population common variants (allele frequency > 1%) in cCREs

ln -s ../reference_files/ALL.GRCh38_sites.20170504.AF0.01.vcf .
ln -s ../reference_files/gencode.v41.annotation.gtf.gz .
ln -s ../compare_TE_nonTE_GWAS/final_info_human_only.all.bed .
ln -s ../human_mouse_cCRE_comparison/final_info_human .
ln -s ../reference_files/hg38_lite.size .

#Find cCREs with protein-coding overlap (for exclusion)
zgrep 'transcript_type "protein_coding"' gencode/gencode.v41.annotation.gtf.gz > gencode_protein_coding_transcripts
awk '{if ($3 == "CDS") print $0}' gencode_protein_coding_transcripts > gencode_protein_coding_transcripts_filtered
cut -f 1,4,5 gencode_protein_coding_transcripts_filtered | sort -k 1,1 -k2,2n | uniq > protein_coding_coords_gencode.bed
bedtools intersect -wao -a final_info_human_only.all.bed -b protein_coding_coords_gencode.bed > overlap_cCREs_protein_coding
awk 'FS="\t" {if ($15 != ".") print $0}' overlap_cCREs_protein_coding > protein_coding_overlap_cCREs
cat ../reference_files/hg38.Ngap.bed protein_coding_coords_gencode.bed > for_exclude_gaps_CDS.bed

#Get common variants in cCREs
perl -lane 'splice @F,10,11;print join("\t",@F)' final_info_human > final_info_human_only.txt
perl -lane '@bed = split /:|-/, $F[0];splice @F, 0,1,@bed;print join("\t",@F) if $F[6] eq "NA"' final_info_human_only.txt > final_info_human_only.noTE.bed
perl -lane '@bed = split /:|-/, $F[0];splice @F, 0,1,@bed;print join("\t",@F) if $F[6] ne "NA" && $. >1' final_info_human_only.txt > final_info_human_only.TE.bed
bedtools intersect -a final_info_human_only.TE.bed -b ALL.GRCh38_sites.20170504.AF0.01.vcf -c > final_info_human_only.TE.countkilovariants.AF0.01.bed
bedtools intersect -a final_info_human_only.noTE.bed -b ALL.GRCh38_sites.20170504.AF0.01.vcf -c > final_info_human_only.noTE.countkilovariants.AF0.01.bed
awk 'OFS="\t" {print $0,$15/($3-$2)*100}' /scratch/xzhuo/for_variation_analysis/final_info_human_only.TE.countkilovariants.AF0.01.bed > observed_varFreq_common_TEs
awk 'OFS="\t" {print $0,$15/($3-$2)*100}' /scratch/xzhuo/for_variation_analysis/final_info_human_only.noTE.countkilovariants.AF0.01.bed > observed_varFreq_common_nonTE
cat observed_varFreq_common_TEs observed_varFreq_common_nonTE > observed_varFreq_common_all
#Remove cCREs that overlap protein-coding sequence
python3 remove_CDS_overlap_cCREs.py observed_varFreq_common_all protein_coding_overlap_cCREs observed_varFreq_common_all_noCDS
python3 remove_CDS_overlap_cCREs.py observed_varFreq_common_TEs protein_coding_overlap_cCREs observed_varFreq_common_TEs_noCDS
python3 remove_CDS_overlap_cCREs.py observed_varFreq_common_nonTE protein_coding_overlap_cCREs observed_varFreq_common_nonTE_noCDS

#Shuffle cCRE positions in the genome while excluding gapped and protein-coding regions
mkdir exclude_CDS_shuffles
mkdir exclude_CDS_overlaps
counter=0
for i in `seq 1 1000`; do
        ((counter+=1))
        bedtools shuffle -i final_info_human_only.all.bed -g hg38_lite.size -noOverlapping -excl for_exclude_gaps_CDS.bed > exclude_CDS_shuffles/shuffle"$i" &
        if [ $counter == 8 ]; then
                wait
                counter=0
        fi
done
#Intersect with human common variants
counter=0
for i in `seq 1 1000`; do
        ((counter+=1))
        bedtools intersect -c -a exclude_CDS_shuffles/shuffle"$i" -b ALL.GRCh38_sites.20170504.AF0.01.vcf | awk 'OFS="\t" {print $0,$15/($3-$2)*100}' > exclude_CDS_overlaps/overlap"$i" &
        if [ $counter == 8 ]; then
                wait
                counter=0
        fi
done
for i in `seq 0 9`; do
        for j in `seq 1 100`; do
                temp=`expr "$i" \* 100`
                counter=`expr "$temp" + "$j"`
                echo exclude_CDS_overlaps/overlap"$counter" >> list_group"$i"
        done
done
#Summarize common variant overlaps
for i in `seq 0 9`; do
        python3 calc_summary_stats_shuffles.py `cat list_group"$i"` summary_stats_shuffle"$i" &
done
wait
cat summary_stats_shuffle0_includeZero > summary_stats_shuffle_all_includeZero
for i in `seq 1 9`; do tail -n+2 summary_stats_shuffle"$i"_includeZero >> summary_stats_shuffle_all_includeZero; done;
cat summary_stats_shuffle0_nonZero > summary_stats_shuffle_all_nonZero
for i in `seq 1 9`; do tail -n+2 summary_stats_shuffle"$i"_nonZero >> summary_stats_shuffle_all_nonZero; done;

#Check flanking regions
bedtools flank -pct -l 1.0 -r 0 -i final_info_human_only.all.bed -g /bar/genomes/hg38/hg38_lite.size > flank_left_cCREs.bed
bedtools flank -pct -r 1.0 -l 0 -i final_info_human_only.all.bed -g /bar/genomes/hg38/hg38_lite.size > flank_right_cCREs.bed
bedtools intersect -c -a flank_left_cCREs.bed -b ALL.GRCh38_sites.20170504.AF0.01.vcf | awk 'OFS="\t" {print $0,$15/($3-$2)*100}' > flank_left_cCREs_overlap
bedtools intersect -c -a flank_right_cCREs.bed -b ALL.GRCh38_sites.20170504.AF0.01.vcf | awk 'OFS="\t" {print $0,$15/($3-$2)*100}' > flank_right_cCREs_overlap
#Remove protein coding overlap
bedtools intersect -a flank_left_cCREs.bed -b ../background_check/protein_coding_coords_gencode.bed > protein_coding_overlap_flank_left
bedtools intersect -a flank_right_cCREs.bed -b ../background_check/protein_coding_coords_gencode.bed > protein_coding_overlap_flank_right
python3 remove_CDS_overlap_cCREs.py flank_left_cCREs_overlap protein_coding_overlap_flank_left flank_left_cCREs_overlap_noCDS
python3 remove_CDS_overlap_cCREs.py flank_right_cCREs_overlap protein_coding_overlap_flank_right flank_right_cCREs_overlap_noCDS

#Get flanking regions of shuffled regions
mkdir exclude_CDS_shuffles_flanks
mkdir exclude_CDS_overlaps_flanks
counter=0
for i in `seq 1 1000`; do
        ((counter+=1))
        bedtools flank -pct -l 1.0 -r 0 -i exclude_CDS_shuffles/shuffle"$i" -g /bar/genomes/hg38/hg38_lite.size > exclude_CDS_shuffles_flanks/left_flank"$i" &
        bedtools flank -pct -r 1.0 -l 0 -i exclude_CDS_shuffles/shuffle"$i" -g /bar/genomes/hg38/hg38_lite.size > exclude_CDS_shuffles_flanks/right_flank"$i" &
        if [ $counter == 4 ]; then
                wait
                counter=0
        fi
done
counter=0
for i in `seq 1 1000`; do
        ((counter+=1))
        bedtools intersect -c -a exclude_CDS_shuffles_flanks/left_flank"$i" -b ALL.GRCh38_sites.20170504.AF0.01.vcf | awk 'OFS="\t" {print $0,$15/($3-$2)*100}' > exclude_CDS_overlaps_flanks/left_overlap"$i" &
        bedtools intersect -c -a exclude_CDS_shuffles_flanks/right_flank"$i" -b ALL.GRCh38_sites.20170504.AF0.01.vcf | awk 'OFS="\t" {print $0,$15/($3-$2)*100}' > exclude_CDS_overlaps_flanks/right_overlap"$i" &
        if [ $counter == 4 ]; then
                wait
                counter=0
        fi
done
#Remove CDS overlapping flanks
counter=0
for i in `seq 1 1000`; do
        ((counter+=1))
        bedtools intersect -a exclude_CDS_shuffles_flanks/left_flank"$i" -b ../background_check/protein_coding_coords_gencode.bed > exclude_CDS_shuffles_flanks/CDS_overlap_left_flank"$i" &
        bedtools intersect -a exclude_CDS_shuffles_flanks/right_flank"$i" -b ../background_check/protein_coding_coords_gencode.bed > exclude_CDS_shuffles_flanks/CDS_overlap_right_flank"$i" &
        if [ $counter == 4 ]; then
                wait
                counter=0
        fi
done
counter=0
for i in `seq 1 1000`; do
        ((counter+=1))
        python3 remove_CDS_overlap_cCREs.py exclude_CDS_overlaps_flanks/left_overlap"$i" exclude_CDS_shuffles_flanks/CDS_overlap_left_flank"$i" exclude_CDS_overlaps_flanks/left_overlap"$i"_noCDS &
        python3 remove_CDS_overlap_cCREs.py exclude_CDS_overlaps_flanks/right_overlap"$i" exclude_CDS_shuffles_flanks/CDS_overlap_right_flank"$i" exclude_CDS_overlaps_flanks/right_overlap"$i"_noCDS &
        if [ $counter == 4 ]; then
                wait
                counter=0
        fi
done
mkdir results_flanks
for dir in left right; do
        rm list_group*
        for i in `seq 0 9`; do
                for j in `seq 1 100`; do
                        temp=`expr "$i" \* 100`
                        counter=`expr "$temp" + "$j"`
                        echo exclude_CDS_overlaps_flanks/"$dir"_overlap"$counter"_noCDS >> list_group"$i"
                done
        done
        for i in `seq 0 9`; do
                python3 calc_summary_stats_shuffles.py `cat list_group"$i"` results_flanks/summary_stats_flank_"$dir""$i" &
        done
        wait
done
for dir in left right; do
        cat results_flanks/summary_stats_flank_"$dir"0_includeZero > results_flanks/summary_stats_flank_"$dir"_all_includeZero
        cat results_flanks/summary_stats_flank_"$dir"0_nonZero > results_flanks/summary_stats_flank_"$dir"_all_nonZero
        for i in `seq 1 9`; do tail -n+2 results_flanks/summary_stats_flank_"$dir""$i"_includeZero >> results_flanks/summary_stats_flank_"$dir"_all_includeZero; done;
        for i in `seq 1 9`; do tail -n+2 results_flanks/summary_stats_flank_"$dir""$i"_nonZero >> results_flanks/summary_stats_flank_"$dir"_all_nonZero; done;
done

#Calculate empirical p-values for cCRE vs. flanks (observed vs. shuffled)
python3 calc_empirical_pvals_flanks_nonZero.py observed_varFreq_common_all_noCDS flank_left_cCREs_overlap_noCDS flank_right_cCREs_overlap_noCDS summary_stats_shuffle_all_nonZero results_flanks/summary_stats_flank_left_all_nonZero results_flanks/summary_stats_flank_right_all_nonZero empirical_pvals_flanks_nonZero
python3 calc_empirical_pvals_flanks_nonZero_v2.py observed_varFreq_common_all_noCDS flank_left_cCREs_overlap_noCDS flank_right_cCREs_overlap_noCDS summary_stats_shuffle_all_nonZero results_flanks/summary_stats_flank_left_all_nonZero results_flanks/summary_stats_flank_right_all_nonZero v2_empirical_pvals_flanks_nonZero

mkdir forR_plotting
python3 consolidate_observed_flanks_forR.py observed_varFreq_common_all_noCDS flank_left_cCREs_overlap_noCDS flank_right_cCREs_overlap_noCDS forR_plotting/flank_varFreq

for i in `seq 1 100`; do cat exclude_CDS_overlaps/overlap"$i" >> shuffle_varFreq_common_100combined_noCDS; done;

