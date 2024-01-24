#Commands used to profile TE-derived cCREs

#Create bed files for each TE subfamily
mkdir separate_TEsubfamily_bed
cd separate_TEsubfamily_bed
perl -lane '{push @{$h{$F[3]}}, $_}END{for $key (keys %h){open my $fh, ">", $key.".bed";for $i (@{$h{$key}}){print $fh "$i"}close()}}' ../reAnno_hg38.bed
cd ..

mkdir forR_plotting
#Overlap with cCREs (all cCREs and K562)
sort -k 1,1 -k2,2n hg38_ccres.bed > hg38_ccres_sorted.bed
bedtools intersect -f 0.5 -wao -sorted -a hg38_ccres_sorted.bed -b reAnno_hg38.bed > allHg38ccre_allTE_overlap
awk '{if ($7!=".") print $0}' allHg38ccre_allTE_overlap > allHg38ccre_allTE_overlap_only
#Add header
sed -i '1i chr\tstart\tend\tnot_sure\tID\tcCRE-type\tTEchr\tTEstart\tTEend\trepName\trepClass/Family\tstrand\talignment\tbpOverlap' allHg38ccre_allTE_overlap
sed -i '1i chr\tstart\tend\tnot_sure\tID\tcCRE-type\tTEchr\tTEstart\tTEend\trepName\trepClass/Family\tstrand\talignment\tbpOverlap' allHg38ccre_allTE_overlap_only

#Subfamily enrichments
python3 get_subfamily_enrichments_cCRE_overlap.py reAnno_hg38.bed allHg38ccre_allTE_overlap subfamily_cCRE_enrichments_all > subfamily_cCRE_overlap_numbers_all
python3 get_subfamily_enrichments_cCRE_overlap_combined.py reAnno_hg38.bed allHg38ccre_allTE_overlap subfamily_cCRE_enrichments_all_combined > subfamily_cCRE_overlap_numbers_all_combined
#Subfamily percentages
awk '{print $1, $2, $4/$2, $6/$2, $8/$2, $10/$2, $12/$2}' <(tail -n+2 subfamily_cCRE_overlap_numbers_all_combined) | sed 's/\s/\t/g' > subfamily_cCRE_percentages_combined
sed -i '1 i\Subfamily\tSubfamily_count\tpELS_percent\tdELS_percent\tPLS_percent\tCTCF-only_percent\tDNase-H3K4me3_percent' subfamily_cCRE_percentages_combined

#Separate cCREs and subfamilies to do fishers tests
#cCRE categories are: pELS, dELS, PLS, DNase-H3K4me3, CTCF-only
mkdir separate_cCRE_files/all_hg38
bedtools sort -i <(grep "pELS" hg38_ccres.bed) > separate_cCRE_files/all_hg38/hg38_ccres_pELS.bed
bedtools sort -i <(grep "dELS" hg38_ccres.bed) > separate_cCRE_files/all_hg38/hg38_ccres_dELS.bed
bedtools sort -i <(grep "PLS" hg38_ccres.bed) > separate_cCRE_files/all_hg38/hg38_ccres_PLS.bed
bedtools sort -i <(grep "DNase-H3K4me3" hg38_ccres.bed) > separate_cCRE_files/all_hg38/hg38_ccres_DNase-H3K4me3.bed
bedtools sort -i <(grep "CTCF-only" hg38_ccres.bed) > separate_cCRE_files/all_hg38/hg38_ccres_CTCF-only.bed
#Perform fishers tests using bedtools
#All hg38 cCREs
mkdir fishers_hg38_cCREs_overlap_TEsubfamily
ls separate_cCRE_files/all_hg38/ | cut -d"_" -f 3 | cut -d"." -f 1 > list_ccre_types
ls separate_TEsubfamily_bed/ | cut -d"." -f 1 > list_TE_subfamilies
mkdir TEsubfamilies_grouped
python3 split_TEsubfamily_list.py list_TE_subfamilies 12 .  TEsubfamilies_grouped/TEsubfamily
for i in `cat list_ccre_types`; do
        for group in `ls TEsubfamilies_grouped/`; do
                for j in `cat TEsubfamilies_grouped/"$group"`; do
                        bedtools fisher -f 0.5 -sorted -a separate_cCRE_files/all_hg38/hg38_ccres_"$i".bed -b separate_TEsubfamily_bed/"$j".bed -g hg38.chrom.sort.sizes > fishers_hg38_cCREs_overlap_TEsubfamily/"$j"_"$i".txt &
                done
                wait
        done
done
python3 summarize_fishers_cCRE_subfamily_overlaps.py fishers_hg38_cCREs_overlap_TEsubfamily/ allHg38ccre_allTE_fishers
#Look for interesting subfamilies
cat <(head -n1 allHg38ccre_allTE_fishers ) <(tail -n+2 allHg38ccre_allTE_fishers | sort -k4 -g -r ) > allHg38ccre_allTE_fishers_sortNum
cat <(head -n1 allHg38ccre_allTE_fishers ) <(awk '{if ($6<0.05) print $0}' allHg38ccre_allTE_fishers_sortNum) > allHg38ccre_allTE_fishers_sortNum_enriched

#Separate TE classes
mkdir separate_TEclasses
cd separate_TEclasses
perl -lane '{@class=split /\//, $F[4];push @{$h{$class[0]}}, $_}END{for $key (keys %h){open my $fh, ">", $key.".bed";for $i (@{$h{$key}}){print $fh "$i"}close()}}' ../reAnno_hg38.bed
cd ..
mkdir separate_TEclasses_bed
sort -k 1,1 -k2,2n <(cat separate_TEclasses/LTR*.bed) > separate_TEclasses_bed/LTR.bed
sort -k 1,1 -k2,2n <(cat separate_TEclasses/DNA*.bed) > separate_TEclasses_bed/DNA.bed
sort -k 1,1 -k2,2n <(cat separate_TEclasses/SINE*.bed) > separate_TEclasses_bed/SINE.bed
sort -k 1,1 -k2,2n <(cat separate_TEclasses/LINE*.bed) > separate_TEclasses_bed/LINE.bed
sort -k 1,1 -k2,2n <(cat separate_TEclasses/Retroposon.bed) > separate_TEclasses_bed/SVA.bed
#Fishers tests
mkdir fishers_hg38_cCREs_overlap_TEclasses
for i in `cat list_ccre_types`; do
        for j in `ls separate_TEclasses_bed/ | cut -d"." -f 1`; do
                bedtools fisher -f 0.5 -sorted -a separate_cCRE_files/all_hg38/hg38_ccres_"$i".bed -b separate_TEclasses_bed/"$j".bed -g hg38.chrom.sort.sizes > fishers_hg38_cCREs_overlap_TEclasses/"$j"_"$i".txt
        done
done
python3 summarize_fishers_cCRE_subfamily_overlaps.py fishers_hg38_cCREs_overlap_TEclasses/ allHg38ccre_allClasses_fishers
#Overlap and enrichment
python3 get_class_enrichments_cCRE_overlap_combined.py reAnno_hg38.bed allHg38ccre_allTE_overlap classLevel_cCRE_enrichments_all_combined > classLevel_cCRE_overlap_numbers_all_combined

#Cell/tissue type overlaps
#Get Cell/tissue types with full cCRE classification
for i in `ls ccre_files_human/`; do
        echo "$i" > temp
        zcat ccre_files_human/"$i" | head -n1 | cut -f 11 >> temp
        cat temp | tr -s '\n' '\t' | awk '{if ($2=="All-data/Full-classification") print $0}' >> list_full_classification_human_ccre_files
done
mkdir full_classification_human_ccre_files
cd full_classification_human_ccre_files/
for i in `cut -f 1 ../list_full_classification_human_ccre_files`; do ln -s ../ccre_files_human/"$i" . ; done;
cd ..
#Get separate cCREs in each cell/tissue type
mkdir full_classification_human_ccre_files_separateCREs
for i in `cut -f 1 list_full_classification_human_ccre_files | sed 's/\.bed\.gz//g'`; do
        bedtools sort -i <(zgrep "pELS" ccre_files_human/"$i".bed.gz) > full_classification_human_ccre_files_separateCREs/"$i"_pELS.bed
        bedtools sort -i <(zgrep "dELS" ccre_files_human/"$i".bed.gz) > full_classification_human_ccre_files_separateCREs/"$i"_dELS.bed
        bedtools sort -i <(zgrep "PLS" ccre_files_human/"$i".bed.gz) > full_classification_human_ccre_files_separateCREs/"$i"_PLS.bed
        bedtools sort -i <(zgrep "DNase-H3K4me3" ccre_files_human/"$i".bed.gz) > full_classification_human_ccre_files_separateCREs/"$i"_DNase-H3K4me3.bed
        bedtools sort -i <(zgrep "CTCF-only" ccre_files_human/"$i".bed.gz) > full_classification_human_ccre_files_separateCREs/"$i"_CTCF-only.bed
        bedtools sort -i <(zgrep "DNase-only" ccre_files_human/"$i".bed.gz) > full_classification_human_ccre_files_separateCREs/"$i"_DNase-only.bed
        bedtools sort -i <(zgrep "CTCF" ccre_files_human/"$i".bed.gz) > full_classification_human_ccre_files_separateCREs/"$i"_CTCF-bound.bed
done

#Get subfamily-TE class association
#Restrict to main 4 classes of DNA, LINE, SINE, LTR
rm temp
for i in `ls separate_TEsubfamily_bed`; do head -n1 separate_TEsubfamily_bed/"$i" >> temp; done;
cut -f 4,5 temp | grep -e "DNA" -e "LTR" -e "LINE" -e "SINE" -e "SVA" > list_TE_classes
grep "SINE" list_TE_classes > list_SINE_subfamilies
grep "DNA" list_TE_classes > list_DNA_subfamilies
grep "LTR" list_TE_classes > list_LTR_subfamilies
grep "LINE" list_TE_classes > list_LINE_subfamilies
grep "SVA" list_TE_classes > list_SVA_subfamilies
cut -f 4,5 temp | sort -g | uniq | grep -v "DNA" | grep -v "LTR" | grep -v "LINE" | grep -v "SINE" | grep -v "SVA" > list_TE_excluded_classes

#Get significant subfamilies from fishers tests for enrichment plotting
tail -n+2 allHg38ccre_allTE_fishers_sortNum | awk '{if ($7<0.05) print $0}' | awk '{if ($3>50) print $0}' | cut -f 1 | sort -g | uniq > sigFishers_subfamilies
python3 create_subfamily_cCRE_enrichment_matrix_forR.py subfamily_cCRE_enrichments_all_combined list_TE_classes sigFishers_subfamilies bases sigFishers_subfamilies_enrichBases_forR
#Get class enrichments for plotting
python3 create_class_cCRE_enrichment_matrix_forR.py classLevel_cCRE_enrichments_all_combined bases classLevel_cCRE_enrichments_all_combined_forR
#Get subfamily enrichment distributions for each main class
python3 create_subfamily_cCRE_enrichment_matrix_forR.py subfamily_cCRE_enrichments_all_combined list_TE_classes list_TE_subfamilies bases subfamily_cCRE_enrichments_all_combined_forR

#Do intersects of TE subfamilies and 25 full classification cell/tissue types
mkdir overlaps_TEsubfamilies_fullClassCCREs
mkdir fishers_TEsubfamilies_fullClassCCREs
for i in `cat list_TE_subfamilies`; do
        mkdir overlaps_TEsubfamilies_fullClassCCREs/"$i"
        mkdir fishers_TEsubfamilies_fullClassCCREs/"$i"
        for j in `ls full_classification_human_ccre_files_separateCREs | sed 's/\.bed//g'`; do
                bedtools intersect -F 0.5 -sorted -wao -a separate_TEsubfamily_bed/"$i".bed -b full_classification_human_ccre_files_separateCREs/"$j".bed > overlaps_TEsubfamilies_fullClassCCREs/"$i"/"$j"."$i" &
                bedtools fisher -F 0.5 -sorted -a separate_TEsubfamily_bed/"$i".bed -b full_classification_human_ccre_files_separateCREs/"$j".bed -g hg38.chrom.sort.sizes > fishers_TEsubfamilies_fullClassCCREs/"$i"/"$j".txt &
        done
        wait
done
#Get subfamily enrichments
mkdir enrichments_TEsubfamilies_fullClassCCREs
mkdir fishers_TEsubfamilies_fullClassCCREs_summaries
for i in `cat list_TE_subfamilies`; do
        python3 get_base_enrichment_from_overlap_single_condition_v2.py overlaps_TEsubfamilies_fullClassCCREs/"$i"/ full_classification_human_ccre_files/ enrichments_TEsubfamilies_fullClassCCREs/"$i"_enrichment_base_full_classification_summary &
        python3 summarize_fishers_subfamily_cCRE_overlaps.py fishers_TEsubfamilies_fullClassCCREs/"$i"/ fishers_TEsubfamilies_fullClassCCREs_summaries/"$i"_fishers_full_classification_summary &
        wait
done
#Combine enrichments
mkdir forR_plotting
python3 combine_enrichments_TEsubfamilies_fullClassCCREs_forR.py enrichments_TEsubfamilies_fullClassCCREs/ 50 forR_plotting/enrichments_TEsubfamilies_fullClassCCREs
python3 combine_enrichments_TEsubfamilies_fullClassCCREs_forR_wClass.py enrichments_TEsubfamilies_fullClassCCREs/ list_TE_classes 0 forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass
#Get distribution of the subfamily enrichments for each TE class
cat <(tail -n+2 forR_plotting/enrichments_TEsubfamilies_fullClassCCREs_wClass_*) > temp
python3 get_distribution_subfamily_celltype_enrichments.py temp 1 forR_plotting/distribution_TEclass_subfamily_celltype_enrichments

#Get element level number of overlaps for all subfamilies
mkdir element_level_overlaps
for group in `ls TEsubfamilies_grouped/`
do
        for i in `cat TEsubfamilies_grouped/"$group"`
        do
                python3 get_ccre_overlaps_per_element_v2.py overlaps_TEsubfamilies_fullClassCCREs/"$i"/ element_level_overlaps/"$i"_element_overlaps_fullClassCCREs &
        done
        wait
done

#Get subfamily overlap number distributions for each main class
#Overlap counts
python3 create_subfamily_cCRE_enrichment_matrix_forR.py subfamily_cCRE_overlap_numbers_all_combined list_TE_classes list_TE_subfamilies count subfamily_cCRE_overlap_numbers_all_combined_forR_count
sed -i 's/FoldEnrich/Number/g' subfamily_cCRE_overlap_numbers_all_combined_forR_count
grep -v "Subfamily_count" subfamily_cCRE_overlap_numbers_all_combined_forR_count > temp && mv temp subfamily_cCRE_overlap_numbers_all_combined_forR_count
#Overlap bases
python3 create_subfamily_cCRE_enrichment_matrix_forR.py subfamily_cCRE_overlap_numbers_all_combined list_TE_classes list_TE_subfamilies bases subfamily_cCRE_overlap_numbers_all_combined_forR_bases
sed 's/FoldEnrich/Number/g' subfamily_cCRE_overlap_numbers_all_combined_forR_bases | grep -v "Subfamily_count" > temp && mv temp subfamily_cCRE_overlap_numbers_all_combined_forR_bases

#Get overlap numbers for each subfamily across 25 full classification cell/tissue types
#Combine overlap numbers
python3 combine_overlaps_TEsubfamilies_fullClassCCREs_forR_wClass.py enrichments_TEsubfamilies_fullClassCCREs/ list_TE_classes 0 forR_plotting/overlapCounts_TEsubfamilies_fullClassCCREs_wClass
#Plot using plotting_subfamily_cCRE_overlap_numbers.R

#Look at fisher's test significant and enriched (by ratio) subfamilies
tail -n+1 allHg38ccre_allTE_fishers_sortNum_enriched | cut -f 1 | sort -g | uniq > sigFishers_enriched_subfamilies
python3 create_subfamily_cCRE_enrichment_matrix_forR.py subfamily_cCRE_enrichments_all_combined list_TE_classes sigFishers_enriched_subfamilies bases sigFishers_enriched_subfamilies_enrichBases_forR

#Get overall proportion of TE-cCRE overlap
python3 calc_cCRE_proportions.py allHg38ccre_allTE_overlap proportion_cCREs_classes
grep -v "All_REs" proportion_cCREs_classes | grep -v "Other" > proportion_cCREs_classes_justTEs
#Get overall TE proportion in the genome
python3 calc_TE_proportions.py reAnno_hg38.bed proportion_genomic_TEs

#Get overall proportions of TE-cCRE overlap for each full classified cell type
mkdir full_classification_human_ccre_files
for i in `ls ../full_classification_human_ccre_files/ | sed 's/.gz//g'`; do zgrep -v "Low-DNase" ../full_classification_human_ccre_files/"$i".gz | sort -k 1,1 -k2,2n > full_classification_human_ccre_files/"$i"; done;
ls full_classification_human_ccre_files/ | sed 's/.bed//g' > list_full_classification_celltypes
mkdir overlaps_allTEs_fullClassCCREs
mkdir proportions_fullClassCCREs
for i in `cat list_full_classification_celltypes`; do
        bedtools intersect -f 0.5 -wao -sorted -a full_classification_human_ccre_files/"$i".bed -b reAnno_hg38.bed > overlaps_allTEs_fullClassCCREs/"$i".cCRE_overlap
        python3 calc_cCRE_proportions_celltype.py overlaps_allTEs_fullClassCCREs/"$i".cCRE_overlap proportions_fullClassCCREs/"$i"_proportion_cCREs
done
python3 split_proportions_fullClassCells_byCCRE.py proportions_fullClassCCREs/ forR_plotting/proportions_TEclass_fullClassCCREs

#Get number of elements per TE subfamily
cut -f 4 reAnno_hg38.bed | sort | uniq -c > number_elements_perSubfamily
echo -e "Subfamily\tNumber" > number_elements_perSubfamily_org
awk 'OFS="\t" {print $2, $1}' number_elements_perSubfamily >> number_elements_perSubfamily_org

#Get numbers for any overlap (expecially bp of overlap)
mkdir any_overlap
bedtools intersect -wao -sorted -a hg38_ccres_sorted.bed -b reAnno_hg38.bed > any_overlap/allHg38ccre_allTE_overlap
awk '{if ($7!=".") print $0}' any_overlap/allHg38ccre_allTE_overlap > any_overlap/allHg38ccre_allTE_overlap_only

#Quantify cell type specificity based on full classification cell/tissue types
for f in full_classification_human_ccre_files_separateCREs/*.bed; do perl -slpe '$_.="\t$f"' -- -f=$f $f >> ENCODTE.25all.bed;done
sort -k 1,1 -k2,2n ENCODTE.25all.bed > ENCODTE.25all.sort.bed
perl -lane '$F[-1] =~/(^\S\S\S+?)_/;$F[-1]=$1;@cres=split /,/, $F[9];$F[4]=$cres[0];splice @F,5,-1;print join("\t",@F)' ENCODTE.25all.sort.bed > ENCODTE.25all.pretty.sort.bed
perl -lane 'print $_ if $F[4] eq "PLS"' ENCODTE.25all.pretty.sort.bed > ENCODTE.25.PLS.pretty.sort.bed
perl -lane 'print $_ if $F[4] eq "pELS"' ENCODTE.25all.pretty.sort.bed > ENCODTE.25.pELS.pretty.sort.bed
perl -lane 'print $_ if $F[4] eq "dELS"' ENCODTE.25all.pretty.sort.bed > ENCODTE.25.dELS.pretty.sort.bed
perl -lane 'print $_ if $F[4] eq "CTCF"' ENCODTE.25all.pretty.sort.bed > ENCODTE.25.CTCF.pretty.sort.bed
perl -lane 'print $_ if $F[4] eq "DNase"' ENCODTE.25all.pretty.sort.bed > ENCODTE.25.DNase.pretty.sort.bed

bedtools groupby -i ENCODTE.25.PLS.pretty.sort.bed -g 1-5 -c 6 -o count,distinct > ENCODTE.25.PLS.pretty.groupby.bed
bedtools groupby -i ENCODTE.25.pELS.pretty.sort.bed -g 1-5 -c 6 -o count,distinct > ENCODTE.25.pELS.pretty.groupby.bed
bedtools groupby -i ENCODTE.25.dELS.pretty.sort.bed -g 1-5 -c 6 -o count,distinct > ENCODTE.25.dELS.pretty.groupby.bed
bedtools groupby -i ENCODTE.25.CTCF.pretty.sort.bed -g 1-5 -c 6 -o count,distinct > ENCODTE.25.CTCF.pretty.groupby.bed
bedtools groupby -i ENCODTE.25.DNase.pretty.sort.bed -g 1-5 -c 6 -o count,distinct > ENCODTE.25.DNase.pretty.groupby.bed

perl -lane 'BEGIN{@cells=("A673","AG04450","B_cell","CD14-positive","GM12878","GM23338","H1_embryonic","HCT116","HeLa-S3","HepG2","IMR-90","K562","MCF-7","MM.1S","OCI-LY7","PC-3","PC-9","Panc1","astrocyte","bipolar","fibroblast","hepatocyte","keratinocyte","myotube","neural");print "chr\tstart\tend\tname\tcCRE\ttotal\t".join("\t",@cells)}{@t=split /,/, @F[6];@mat=map {$_~~@t?1:0} @cells; splice @F, 6, 2, @mat;print join("\t",@F)}' ENCODTE.25.PLS.pretty.groupby.bed > ENCODTE.25.PLS.upset.txt
perl -lane 'BEGIN{@cells=("A673","AG04450","B_cell","CD14-positive","GM12878","GM23338","H1_embryonic","HCT116","HeLa-S3","HepG2","IMR-90","K562","MCF-7","MM.1S","OCI-LY7","PC-3","PC-9","Panc1","astrocyte","bipolar","fibroblast","hepatocyte","keratinocyte","myotube","neural");print "chr\tstart\tend\tname\tcCRE\ttotal\t".join("\t",@cells)}{@t=split /,/, @F[6];@mat=map {$_~~@t?1:0} @cells; splice @F, 6, 2, @mat;print join("\t",@F)}' ENCODTE.25.pELS.pretty.groupby.bed > ENCODTE.25.pELS.upset.txt
perl -lane 'BEGIN{@cells=("A673","AG04450","B_cell","CD14-positive","GM12878","GM23338","H1_embryonic","HCT116","HeLa-S3","HepG2","IMR-90","K562","MCF-7","MM.1S","OCI-LY7","PC-3","PC-9","Panc1","astrocyte","bipolar","fibroblast","hepatocyte","keratinocyte","myotube","neural");print "chr\tstart\tend\tname\tcCRE\ttotal\t".join("\t",@cells)}{@t=split /,/, @F[6];@mat=map {$_~~@t?1:0} @cells; splice @F, 6, 2, @mat;print join("\t",@F)}' ENCODTE.25.dELS.pretty.groupby.bed > ENCODTE.25.dELS.upset.txt
perl -lane 'BEGIN{@cells=("A673","AG04450","B_cell","CD14-positive","GM12878","GM23338","H1_embryonic","HCT116","HeLa-S3","HepG2","IMR-90","K562","MCF-7","MM.1S","OCI-LY7","PC-3","PC-9","Panc1","astrocyte","bipolar","fibroblast","hepatocyte","keratinocyte","myotube","neural");print "chr\tstart\tend\tname\tcCRE\ttotal\t".join("\t",@cells)}{@t=split /,/, @F[6];@mat=map {$_~~@t?1:0} @cells; splice @F, 6, 2, @mat;print join("\t",@F)}' ENCODTE.25.CTCF.pretty.groupby.bed > ENCODTE.25.CTCF.upset.txt
perl -lane 'BEGIN{@cells=("A673","AG04450","B_cell","CD14-positive","GM12878","GM23338","H1_embryonic","HCT116","HeLa-S3","HepG2","IMR-90","K562","MCF-7","MM.1S","OCI-LY7","PC-3","PC-9","Panc1","astrocyte","bipolar","fibroblast","hepatocyte","keratinocyte","myotube","neural");print "chr\tstart\tend\tname\tcCRE\ttotal\t".join("\t",@cells)}{@t=split /,/, @F[6];@mat=map {$_~~@t?1:0} @cells; splice @F, 6, 2, @mat;print join("\t",@F)}' ENCODTE.25.DNase.pretty.groupby.bed > ENCODTE.25.DNase.upset.txt


bedtools intersect -a ENCODTE.25.PLS.upset.txt -b reAnno_hg38.bed -loj |perl -lane 'if ($F[-3] eq "."){$F[-3] = "None"}else{@class=split /\//, $F[-3];$F[-3]=$class[0]};print join("\t",@F)' > ENCODTE.25.PLS.TE.txt
bedtools intersect -a ENCODTE.25.pELS.upset.txt -b reAnno_hg38.bed -loj |perl -lane 'if ($F[-3] eq "."){$F[-3] = "None"}else{@class=split /\//, $F[-3];$F[-3]=$class[0]};print join("\t",@F)' > ENCODTE.25.pELS.TE.txt
bedtools intersect -a ENCODTE.25.dELS.upset.txt -b reAnno_hg38.bed -loj |perl -lane 'if ($F[-3] eq "."){$F[-3] = "None"}else{@class=split /\//, $F[-3];$F[-3]=$class[0]};print join("\t",@F)' > ENCODTE.25.dELS.TE.txt
bedtools intersect -a ENCODTE.25.CTCF.upset.txt -b reAnno_hg38.bed -loj |perl -lane 'if ($F[-3] eq "."){$F[-3] = "None"}else{@class=split /\//, $F[-3];$F[-3]=$class[0]};print join("\t",@F)' > ENCODTE.25.CTCF.TE.txt
bedtools intersect -a ENCODTE.25.DNase.upset.txt -b reAnno_hg38.bed -loj |perl -lane 'if ($F[-3] eq "."){$F[-3] = "None"}else{@class=split /\//, $F[-3];$F[-3]=$class[0]};print join("\t",@F)' > ENCODTE.25.DNase.TE.txt

perl -lane '@classes=("class","LTR","LINE","SINE","DNA","Retroposon");$F[-3] = "None" unless $F[-3] ~~ @classes;print join("\t",@F)' ENCODTE.25.PLS.TE.txt > ENCODTE.25.PLS.TE.5classes.txt
perl -lane '@classes=("class","LTR","LINE","SINE","DNA","Retroposon");$F[-3] = "None" unless $F[-3] ~~ @classes;print join("\t",@F)' ENCODTE.25.pELS.TE.txt > ENCODTE.25.pELS.TE.5classes.txt
perl -lane '@classes=("class","LTR","LINE","SINE","DNA","Retroposon");$F[-3] = "None" unless $F[-3] ~~ @classes;print join("\t",@F)' ENCODTE.25.dELS.TE.txt > ENCODTE.25.dELS.TE.5classes.txt
perl -lane '@classes=("class","LTR","LINE","SINE","DNA","Retroposon");$F[-3] = "None" unless $F[-3] ~~ @classes;print join("\t",@F)' ENCODTE.25.CTCF.TE.txt > ENCODTE.25.CTCF.TE.5classes.txt
perl -lane '@classes=("class","LTR","LINE","SINE","DNA","Retroposon");$F[-3] = "None" unless $F[-3] ~~ @classes;print join("\t",@F)' ENCODTE.25.DNase.TE.txt > ENCODTE.25.DNase.TE.5classes.txt

