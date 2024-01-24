#Commands used to investigate GWAS variant enrichment within in TE-derived cCREs

#GWAS files from EBI (https://www.ebi.ac.uk/gwas/docs/file-downloads)
#Date accessed: 07/07/2023
gwas_catalog_trait-mappings_r2023-07-05.tsv
gwas_catalog_v1.0.2-associations_e109_r2023-07-05.tsv
#These were added/last updated April 24, 2018

#dbSNP153 from UCSC
#Date accessed: 07/07/2023
wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb
#Last modified 2019-11-24 (november 24, 2019)

ln -s ../reference_files/hg38.chrom.sizes .

#Get unique parent GWAS terms
cut -f 4-5 gwas_catalog_trait-mappings_r2023-07-05.tsv | tail -n+2 | sort | uniq > list_parent_terms_unique

#Get GWAS variant coordinates if not provided and confirm those that are provided
cut -f 21 gwas_catalog_v1.0.2-associations_e109_r2023-07-05.tsv | grep "rs" | cut -d"-" -f 1 | sort | uniq > all_rs_gwas_snps
cut -f 23,24 gwas_catalog_v1.0.2-associations_e109_r2023-07-05.tsv | awk '{if ($1 == 1 && $2 != "") print "rs"$2}' | sort | uniq >> all_rs_gwas_snps
sort all_rs_gwas_snps | uniq > temp && mv temp all_rs_gwas_snps
bigBedNamedItems -nameFile dbSnp153.bb all_rs_gwas_snps snp_info_gwas
cut -f 2,8,12,13,21,23,24,35,36 gwas_catalog_v1.0.2-associations_e109_r2023-07-05.tsv > cut_gwas_associations
python3 match_coords_snpID.py cut_gwas_associations snp_info_gwas hg38.chrom.sizes filtered_gwas_associations > problems_gwas_entries 2> problems_matching_gwas_dbsnp
sort -k 1,1 -k2,2n filtered_gwas_associations > filtered_gwas_associations_sorted

#Overlap all cCREs with GWAS variants
tail -n+2 ../human_mouse_cCRE_comparison/final_info_human | sed -e 's/:/\t/' -e 's/-/\t/' | cut -f 1-12,24,25 > final_info_human_only.all.bed
bedtools intersect -sorted -wao -a filtered_gwas_associations_sorted -b final_info_human_only.all.bed > overlap_gwas_observed_cCREs
python3 summarize_gwas_ccre_overlap.py overlap_gwas_observed_cCREs gwas_catalog_trait-mappings_r2023-07-05.tsv counts_observed_gwas_cCRE_overlap > errors_summary

#Shuffle cCRE coordinates for permutation testing of random chance to at least as many GWAS variants as observed
ln -s ../reference_files/hg38_lite.size .
ln -s ../reference_files/hg38.Ngap.bed .
mkdir shuffled_cCRE_locations
mkdir shuffled_overlaps
mkdir shuffled_counts
for i in `seq 1 1000`; do
        bedtools shuffle -i final_info_human_only.all.bed -g hg38_lite.size -noOverlapping -excl hg38.Ngap.bed | sort -k 1,1 -k2,2n > shuffled_cCRE_locations/shuffled_cCREs"$i".bed
        bedtools intersect -sorted -wao -a filtered_gwas_associations_sorted -b shuffled_cCRE_locations/shuffled_cCREs"$i".bed > shuffled_overlaps/overlap_gwas_shuffled_cCREs"$i"
        python3 summarize_gwas_ccre_overlap_shuffled.py shuffled_overlaps/overlap_gwas_shuffled_cCREs"$i" gwas_catalog_trait-mappings_r2023-07-05.tsv shuffled_counts/counts_shuffled_gwas_cCRE_overlap"$i"
done
#Calculate empirical p-values based on permutations
python3 calc_enrichment_observed_vs_shuffled.py counts_observed_gwas_cCRE_overlap shuffled_counts/ enrichments_gwas_cCRE_overlap > pval_observed

#Do an additional 100 shuffled coordinates for negative controls (random genomic regions)
mkdir negative_control
for i in `seq 1 100`; do
        bedtools shuffle -i final_info_human_only.all.bed -g hg38_lite.size -noOverlapping -excl hg38.Ngap.bed | sort -k 1,1 -k2,2n > negative_control/negative_shuffled_ccres"$i"
        bedtools intersect -sorted -wao -a filtered_gwas_associations_sorted -b negative_control/negative_shuffled_ccres"$i" > negative_control/negative_shuffled_overlap"$i"
        python3 summarize_gwas_ccre_overlap_shuffled.py negative_control/negative_shuffled_overlap"$i" gwas_catalog_trait-mappings_r2023-07-05.tsv negative_control/negative_shuffled_counts"$i"
done
python3 average_negative_counts.py negative_control/negative_shuffled_counts* negative_shuffled_counts
python3 calc_enrichment_observed_vs_shuffled_negative.py negative_shuffled_counts shuffled_counts/ negative_shuffled_enrichments > pval_negative_shuffled

