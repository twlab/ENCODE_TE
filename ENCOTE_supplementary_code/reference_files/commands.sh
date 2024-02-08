#Commands to generate reference files

#hg38 repeatmasker file
zcat hg38.fa.out.gz | sed '1,3d' | awk '{OFS="\t"; if (!($11 ~ /Simple_repeat/ || $11 ~ /Low_complexity/ || $11 ~ /Satellite/)) print}' | egrep -v '_alt|_random|chrUn' | pigz -p 16 -c > hg38.filt.fa.out.gz
python controlMerge.v3.py <(zcat hg38.filt.fa.out.gz ) tmp2 tmp
awk '{OFS="\t"; print $1, $2 - 1, $3, $4, $5, $6, $7, $8}' tmp2 | pigz -c -p 16 > reAnno_hg38_ENCOTE.txt.gz
zcat reAnno_hg38_ENCOTE.txt.gz | cut -f 1-7 | sort -k 1,1 -k2,2n > reAnno_hg38.bed
#mm10 repeatmasker file
zcat mm10.fa.out.gz | sed '1,3d' | awk '{OFS="\t"; if (!($11 ~ /Simple_repeat/ || $11 ~ /Low_complexity/ || $11 ~ /Satellite/ || $11 ~ /Unknown/)) print}' | egrep -v '_alt|_random|chrUn' | pigz -p 16 -c > mm10.filt.fa.out.gz
python controlMerge.v3.py <(zcat mm10.filt.fa.out.gz ) tmp2 tmp
awk '{OFS="\t"; print $1, $2 - 1, $3, $4, $5, $6, $7, $8}' tmp2 | pigz -c -p 16 > reAnno_mm10_ENCOTE.txt.gz
zcat reAnno_mm10_ENCOTE.txt.gz | cut -f 1-7 | sort -k 1,1 -k2,2n > reAnno_mm10_ENCOTE.bed

#hg38 cCRE file (cell agnostic)
zcat ENCFF924IMH.bed.gz | sort -k 1,1 -k2,2n > hg38_ccres_sorted.bed
#mm10 cCRE file (cell agnostic)
zcat ENCFF904ZZH.bed.gz | sort -k 1,1 -k2,2n > ccres_mm10_sorted.bed

#List of hg38 and mm10 TE subfamilies
cut -f 4 reAnno_hg38.bed | sort | uniq > list_hg38_TE_subfamilies
cut -f 4 reAnno_mm10_ENCOTE.bed | sort | uniq > list_mm10_TE_subfamilies
#List of hg38 TE subfamilies with TE class
cut -f 4,5 reAnno_hg38.bed | grep -e "DNA" -e "LTR" -e "LINE" -e "SINE" -e "SVA" | sort | uniq > list_TE_classes_human

#Get Kimura % divergence (hg38)
perl calcDivergenceFromAlign.pl -s kimura_distance_hg38.txt hg38.fa.align
#Take the 3' end part of LINEs (fragmented in the RepeatMasker align and subsequent divergence file)
python3 convert_rmsk_kimura_distance_output.py kimura_distance_hg38.txt list_hg38_TE_subfamilies hg38_LINE.json kimura_distance_hg38.tsv > subfamily_name_connections_hg38 2> noKimura_distance_TEsubfamilies
#Get Kimura % divergence (mm10)
perl calcDivergenceFromAlign.pl -s kimura_distance_mm10.txt mm10.fa.align
python3 convert_rmsk_kimura_distance_output.py kimura_distance_mm10.txt list_mm10_TE_subfamilies hg38_LINE.json kimura_distance_mm10.tsv > subfamily_name_connections_mm10 2> failed_kimura_div_mm10

#Use Dfam TE annotations
wget https://www.dfam.org/releases/Dfam_3.6/families/Dfam_curatedonly.embl.gz
gunzip Dfam_curatedonly.embl.gz
#Convert from embl to fasta format
tail -n+32 Dfam_curatedonly.embl > dfam.embl
seqret -sequence dfam.embl -sformat1 embl -outseq Dfam_curatedonly.fa -osformat2 pearson

#Get TE subfamily Clades from Dfam
python3 get_clades_dfam.py dfam.embl clades_dfam_all_TEs
#Manually change some of the clade names
sed -i 's/Theria <mammals>/Theria/g' clades_dfam_all_TEs
sed -i 's/Vertebrata <vertebrates>/Vertebrata/g' clades_dfam_all_TEs
sed -i 's/Mus <genus>/Mus/g' clades_dfam_all_TEs

#Get human (hg38) TE subfamily clades, separated by primate/rodent specific and shared
python3 separate_TEs_cladeBased.py list_hg38_TE_subfamilies subfamily_name_connections_hg38 list_primate_human_lineage_clades ordered_clades_human_mouse_lineages clades_dfam_all_TEs primate_human_lineage_TEs
python3 separate_TEs_cladeBased.py list_hg38_TE_subfamilies subfamily_name_connections_hg38 list_shared_human_rodent_clades ordered_clades_human_mouse_lineages clades_dfam_all_TEs shared_human_rodent_TEs_hg38
#Also get mouse (mm10) TE subfamily clades
python3 separate_TEs_cladeBased.py list_mm10_TE_subfamilies subfamily_name_connections_mm10 list_rodent_mouse_lineage_clades ordered_clades_human_mouse_lineages clades_dfam_all_TEs rodent_mouse_lineage_TEs
python3 separate_TEs_cladeBased.py list_mm10_TE_subfamilies subfamily_name_connections_mm10 list_shared_human_rodent_clades ordered_clades_human_mouse_lineages clades_dfam_all_TEs shared_human_rodent_TEs_mm10

#Get subfamily age (kimura divergence) distributions for groups and clades
python3 assign_age.py primate_human_lineage_TEs kimura_distance_hg38.tsv primate_human_lineage_TEs_wAges
python3 assign_age.py shared_human_rodent_TEs_hg38 kimura_distance_hg38.tsv shared_human_rodent_TEs_hg38_wAges
python3 assign_age.py rodent_mouse_lineage_TEs kimura_distance_mm10.tsv rodent_mouse_lineage_TEs_wAges
python3 assign_age.py shared_human_rodent_TEs_mm10 kimura_distance_mm10.tsv shared_human_rodent_TEs_mm10_wAges
#Add TE class and simultaneously remove non-TE repetitive subfamilies
python3 add_TE_classes.py primate_human_lineage_TEs_wAges list_TE_classes_human primate_human_lineage_TEs_wAges_wClass
python3 add_TE_classes.py shared_human_rodent_TEs_hg38_wAges list_TE_classes_human shared_human_rodent_TEs_hg38_wAges_wClass

#TE consensus sequences
#Get RepeatLibrary 20170127 file (RMRBSeqs.fa)
python3 connect_rmsk_align_outfile_repNames.py hg38.fa.out hg38.fa.align rmsk_out_align_repeatName_connections
grep -v "(" rmsk_out_align_repeatName_connections | grep -v "L1" | grep -v "L2" > rmsk_out_align_repeatName_connections_noSimpleRepeat_noLINE
python3 get_consensus_sequences.py list_hg38_TE_subfamilies rmsk_out_align_repeatName_connections hg38TE.consensus.fa RMRBSeqs.fa final_hg38_consensus_seq.fa > missing_hg38_consensus_seq.txt

#Alignments of each TE copy to its subfamily consensus
grep ">" final_hg38_consensus_seq.fa | sed 's/>//g' > list_TE_subfamilies_wConsensus
python3 separate_consensus_sequences.py final_hg38_consensus_seq.fa TEsubfamily_consensus_fastas
mkdir needle_aligns
mkdir alignFastas
mkdir element_to_consensus_align_metrics
mkdir transition_transversion_rates_TEsubfamilies
for i in `cat list_TE_subfamilies_wConsensus`; do
        needle -asequence TEsubfamily_consensus_fastas/"$i"_consensus.fa -bsequence sequences_TEsubfamilies_forNeedle/"$i".fa -gapopen 10 -gapextend 0.5 -outfile needle_aligns/"$i".needle
        mkdir temp_dir
        python3 divide_needleAlign_individualFiles.py needle_aligns/"$i".needle temp_dir/align
        for j in `ls temp_dir`; do
                python3 parse_needleOutput_fasta.py temp_dir/"$j" >> alignFastas/"$i".alignFasta
        done
        rm -r temp_dir
        sed -i 's/_/:/g' alignFastas/"$i".alignFasta
done

#Human population common variants (allele frequency > 1%)
#VCF files downloaded from 1000 genomes project
#ALL.chr*_GRCh38_sites.20170504.vcf.gz files from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/
zcat ALL.*_GRCh38_sites.20170504.vcf.gz|perl -lane 'if(/^#/){print $_}else{%h=split(/=|;/,$F[7]);print "chr".$_ if $h{"AF"}>0.01}' > ALL.GRCh38_sites.20170504.AF0.01.vcf


