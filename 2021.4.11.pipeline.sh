#!/bin/bash
#coding:utf-8
# set several traps
#set -e 
#脚本只要发生错误，就终止执行。
#set -u 
#执行脚本的时候，如果遇到不存在的变量，Bash 默认忽略它。
#set -o pipefail 
#set -e有一个例外情况，就是不适用于管道命令。只要一个子命令失败，整个管道命令就失败，脚本就会终止执行。

#####################################################################################
# Junbo Yang (yang_junbo_hi@126.com)
# Last modified: 7.15 12:34:47 2020 
# Change log : #
# Change log : #
# Usage: sh script.sh path_to_sample

#####################################################################################
# Locations require to change (based to different computers)

# assign reference genome of the species
THREADS="20"
#REF="/media/ruizhi/database/homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa"
#HIAST_INDEX="/media/ruizhi/database/homo_sapiens/hisat2/GRCh38"
#BWA_INDEX="/media/ruizhi/database/homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa"
Bowtie_index="/media/ruizhi/database/mm/bowtie2_grcm38/GRCm38.rRNA"
#####################################################################################


#####################################################################################
#####################################################################################

# change illumina data format

#####################################################################################
mkdir -p Results/1.Quality_control/raw_QC/
mkdir -p Results/1.Quality_control/raw_QC/
# quality control
ls *_R1.fq.gz | while read id; do(mkdir -p Results/1.Quality_control/raw_QC/sfastqc_results/$(basename $id '.fq.gz');fastqc -o Results/1.Quality_control/raw_QC/sfastqc_results/$(basename $id '.fq.gz') $id);done
ls *_R2.fq.gz | while read id; do(mkdir -p Results/1.Quality_control/raw_QC/sfastqc_results/$(basename $id '.fq.gz');fastqc -o Results/1.Quality_control/raw_QC/sfastqc_results/$(basename $id '.fq.gz') $id);done
multiqc Results/1.Quality_control/raw_QC/sfastqc_results/ -o Results/1.Quality_control/multiqc_results
#####################################################################################
#host genome
ls *_R1.fq.gz | while read id; do(i=${id%_*};bowtie2 -p 24 -x ${Bowtie_index} -1 $id -2 ${i}_R2.fq.gz -S $(basename $id '_R1.fq.gz').sam);done
#ls */*_R1.fq.gz | while read id; do(i=${id%_*};bwa mem -t 2 -R '@RG\tID:lib_lane\tPL:Illumina\tLB:lib\tSM:L' -M -k 30 ${BWA_INDEX} $id ${i}_R2.fq.gz -S $(basename $id '_R1.fq.gz').sam);done
#ls *_R1.fq.gz | while read id; do(hisat2 -p ${THREADS} --pen-noncansplice 1000000 -x ${HIAST_INDEX} -1 $id -2 $(basename $id '_R1.fq.gz')_R2.fq.gz -S $(basename $id '_R1.fq.gz').sam); done
#extract unmapped reads
ls *.sam | while read id;do(samtools view -bS -@ ${THREADS} -T ${REF} -h -u -f 12 -F 256 $id > $(basename $id 'sam')unmapped.bam);done
#sort by positon
ls *.unmapped.bam | while read id;do(samtools sort -@ ${THREADS} -n $id $(basename $id 'bam')sorted);done
#bam2fq
ls *.unmapped.sorted.bam | while read id;do(bamToFastq -i $id -fq $(basename $id '.sorted.bam').R1.fq -fq2 $(basename $id '.sorted.bam').R2.fq);done
#assembly
#ls *R1.fq | while read id;do(mkdir -p Results/spades/$(basename $id ".R1.fq");spades.py -k 21,33,55,77 --only-assembler --pe1-1 $id --pe1-2 $(basename $id "R1.fq")R2.fq -o Results/spades/$(basename $id ".unmapped.R1.fq"));done
#filter scaffolds (gt 500)
#ls *R1.fq.gz | while read id;do(python /media/ruizhi/software/min_500.py -i Results/spades/$(basename $id "_R1.fq.gz")/scaffolds.fasta -o Results/spades/$(basename $id "_R1.fq.gz")/scaffolds.500.fasta);done
ls *.unmapped.R1.fq | while read id;do(mkdir -p Results/2.Raw_assembly/;megahit -1 $id -2 $(basename $id '.unmapped.R1.fq').unmapped.R2.fq -o Results/2.Raw_assembly/$(basename $id '.unmapped.R1.fq') --out-prefix $(basename $id '.unmapped.R1.fq'));done

#filter scaffolds (gt 500)
ls *R1.fq | while read id;do(python /media/ruizhi/software/min_500.py -i Results/2.Raw_assembly/$(basename $id '.unmapped.R1.fq')/$(basename $id '.unmapped.R1.fq').contigs.fa -o Results/2.Raw_assembly/$(basename $id ".unmapped.R1.fq")/$(basename $id '.unmapped.R1.fq').contigs.500.fa);done

####################################################################################################
#evaluate
ls *R1.fq | while read id;do(mkdir -p Results/3.Evaluation/$(basename $id '.unmapped.R1.fq');quast.py -o Results/3.Evaluation/$(basename $id '.unmapped.R1.fq')  -t 20 Results/2.Raw_assembly/$(basename $id '.unmapped.R1.fq')/$(basename $id '.unmapped.R1.fq').contigs.fa);done
ls Results/3.Evaluation/*/report.tsv | while read id;do(i=$id.ft;cut -f 2 $id > $i);done
ls Results/3.Evaluation/*/report.tsv | while read id;do(cut -f 1 $id > Results/3.Evaluation/total.head);done
paste -d '\t' Results/3.Evaluation/*/report.tsv.ft > Results/3.Evaluation/total.tsv
paste -d '\t' Results/3.Evaluation/total.head Results/3.Evaluation/total.tsv > Results/3.Evaluation/total.report.xls
sed -i 's/.contigs//g' Results/3.Evaluation/total.report.xls
rm  Results/3.Evaluation/*/report.tsv.ft
###################################################################################################################
#orf prediction
#ls *R1.fq.gz | while read id;do(mkdir -p Results/orf_prediction/$(basename $id "_R1.fq.gz");prodigal -p meta -a Results/orf_prediction/$(basename $id "_R1.fq.gz")/$(basename $id "_R1.fq.gz").protein_seq.fasta -d Results/orf_prediction/$(basename $id "_R1.fq.gz")/$(basename $id "_R1.fq.gz").nucleotide_seq.fasta -o Results/orf_prediction/$(basename $id "_R1.fq.gz")/$(basename $id "_R1.fq.gz").gbk -s Results/orf_prediction/$(basename $id "_R1.fq.gz")/$(basename $id "_R1.fq.gz").poteintial.stat -i Results/spades/$(basename $id "_R1.fq.gz")/scaffolds.500.fasta);done
ls *R1.fq | while read id;do(mkdir -p Results/6.Function_annotation/;prokka --cpus 10 Results/2.Raw_assembly/$(basename $id '.unmapped.R1.fq')/$(basename $id '.unmapped.R1.fq').contigs.500.fa --outdir Results/4.Raw_annotations/$(basename $id '.unmapped.R1.fq')/ --prefix $(basename $id ".unmapped.R1.fq") --metagenome --kingdom Bacteria);done
#non_redundunt_gene
mkdir -p Results/5.Non_redundunt_geneset
cat Results/4.Raw_annotations/*/*.faa > Results/5.Non_redundunt_geneset/total.redundunt.faa.fa
cat Results/4.Raw_annotations/*/*.ffn > Results/5.Non_redundunt_geneset/total.redundunt.ffn.fa
#cd-hit -M 8000 -i Results/5.Non_redundunt_geneset/total.redundunt.faa -o Results/5.Non_redundunt_geneset/non_redundunt.faa.fa
cd-hit -M 8000 -T 30 -i Results/5.Non_redundunt_geneset/total.redundunt.faa.fa -o Results/5.Non_redundunt_geneset/non_redundunt.nt.fa -c 0.9 -aS 0.9 -d 0
python /media/ruizhi/software/seq_format.py -i Results/5.Non_redundunt_geneset/total.redundunt.ffn.fa -o Results/5.Non_redundunt_geneset/total.redundunt.format.ffn.fa
python /media/ruizhi/software/seq_format.py -i Results/5.Non_redundunt_geneset/non_redundunt.nt.fa -o Results/5.Non_redundunt_geneset/non_redundunt.nt.format.fa

python /media/ruizhi/software/ffn.non.py -p Results/5.Non_redundunt_geneset/non_redundunt.nt.format.fa -f Results/5.Non_redundunt_geneset/total.redundunt.format.ffn.fa -o Results/5.Non_redundunt_geneset/non_redundunt.fnn.format.fa
#############################################################################################################

cd Results/

#CARD
mkdir -p 6.Function_annotation/CARD
ls 5.Non_redundunt_geneset/non_redundunt.nt.format.fa | while read id;do(echo -e "Query_id\tSubject_id\tIdentity\tAlign_length\tMismatch\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE-value\tScore" > head.xls; diamond blastp -p 10 -d /media/ruizhi/database/CARD/CARD -e 1e-5 -k 1 --sensitive -q $id -o 6.Function_annotation/CARD/$(basename $id).card.m8; cat 6.Function_annotation/CARD/head.xls 6.Function_annotation/CARD/$(basename $id).card.m8 > 6.Function_annotation/CARD/$(basename $id).card.xls);done
cd 6.Function_annotation/CARD
python /media/ruizhi/database/CARD/script/card_annotation.py -d /media/ruizhi/database/CARD/script/aro.drug.dict -a /media/ruizhi/database/CARD/script/aro.dict -r /media/ruizhi/database/CARD/script/aro.Re.dict -b non_redundunt.nt.format.fa.card.m8 -o non_redundunt.nt.CARD.xls
rm head.xls
cd ../../
#COG
mkdir -p 6.Function_annotation/COG
ls 5.Non_redundunt_geneset/non_redundunt.nt.format.fa | while read id;do(echo -e "Query_id\tSubject_id\tIdentity\tAlign_length\tMismatch\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE-value\tScore" > 6.Function_annotation/COG/head.xls; diamond blastp -p 10 -d /media/ruizhi/database/COG/cog -e 1e-5 -k 1 --sensitive -q $id -o 6.Function_annotation/COG/$(basename $id).cog.m8; cat 6.Function_annotation/COG/head.xls 6.Function_annotation/COG/$(basename $id).cog.m8 > 6.Function_annotation/COG/$(basename $id).cog.xls);done
cd 6.Function_annotation/COG
python /media/ruizhi/database/COG/script/cog_annotation.py -d /media/ruizhi/database/COG/whog.dict -f /media/ruizhi/database/COG/fun.dict -b non_redundunt.nt.format.fa.cog.m8 -o  non_redundunt.nt.COG.annotation
rm head.xls
cd ../../
#PHI
mkdir -p 6.Function_annotation/PHI
ls 5.Non_redundunt_geneset/non_redundunt.nt.format.fa | while read id;do(echo -e "Query_id\tSubject_id\tIdentity\tAlign_length\tMismatch\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE-value\tScore" > 6.Function_annotation/PHI/head.xls; diamond blastp -p 20 -d /media/ruizhi/database/PHI/phi.dmnd -e 1e-5 -k 1 --sensitive -q $id -o 6.Function_annotation/PHI/$(basename $id).phi.m8; cat 6.Function_annotation/PHI/head.xls 6.Function_annotation/PHI/$(basename $id).phi.m8 > 6.Function_annotation/PHI/$(basename $id).phi.xls);done
rm 6.Function_annotation/PHI/head.xls
#prophase
mkdir -p 6.Function_annotation/prophase
ls 5.Non_redundunt_geneset/non_redundunt.nt.format.fa | while read id;do(echo -e "Query_id\tSubject_id\tIdentity\tAlign_length\tMismatch\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE-value\tScore" > 6.Function_annotation/prophase/head.xls; diamond blastp -p 10 -d /media/ruizhi/database/phase/pro -e 1e-5 -k 1 --sensitive -q $id -o 6.Function_annotation/prophase/$(basename $id).phase.m8; cat 6.Function_annotation/prophase/head.xls 6.Function_annotation/prophase/$(basename $id).phase.m8 > 6.Function_annotation/prophase/$(basename $id).phase.xls);done
python /media/ruizhi/database/phase/pro_annotation.py -d /media/ruizhi/database/phase/prophase.dict -b 6.Function_annotation/prophase/non_redundunt.nt.format.fa.phase.m8 -o 6.Function_annotation/prophase/non_redundunt.nt.prophase.xls
rm 6.Function_annotation/prophase/head.xls
#swissprot
mkdir 6.Function_annotation/swissprot
ls 5.Non_redundunt_geneset/non_redundunt.nt.format.fa | while read id;do(echo -e "Query_id\tSubject_id\tIdentity\tAlign_length\tMismatch\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE-value\tScore" > 6.Function_annotation/swissprot/head.xls; diamond blastp -p 10 -d /media/ruizhi/database/swissprot/swissprot -e 1e-5 -k 1 --sensitive -q $id -o 6.Function_annotation/swissprot/$(basename $id).swissprot.m8; cat 6.Function_annotation/swissprot/head.xls 6.Function_annotation/swissprot/$(basename $id).swissprot.m8 > 6.Function_annotation/swissprot/$(basename $id).swissprot.xls);done
rm 6.Function_annotation/swissprot/head.xls
#NR
mkdir 6.Function_annotation/NR/
split -l 200000 -a 3 -d 5.Non_redundunt_geneset/non_redundunt.fnn.format.fa  6.Function_annotation/NR/non_redundunt.ffn.format_chunk_
cd 6.Function_annotation/NR/
ls *chunk_0* | grep -v m8 | while read id;do(blastn -db /media/ruizhi/database/NR.database/nt_format -outfmt 6 -num_threads 30 -max_target_seqs 1 -query $id -out $(basename $id).NR.m8);done
cat *.m8 > NR.annotation
cd ../../
#GO

#CAZy
# python /media/ruizhi/database/CAZy/tools/run_dbcan/run_dbcan.py non_redundunt.nt.fa protein --out_pre non_redundunt --out_dir non_redundunt
#  python /media/ruizhi/software/stast.py -i non_redundunt.CAZy.out.xls -o stast.xls
#eggNOG
mkdir 6.Function_annotation/eggNOG
split -l 200000 -a 3 -d 5.Non_redundunt_geneset/non_redundunt.nt.format.fa 6.Function_annotation/eggNOG/non_redundunt.nt.format.chunk_
ls 6.Function_annotation/eggNOG/*.chunk_0* | while read id;do(/media/ruizhi/software/eggnog-mapper/emapper.py -m diamond --no_annot --no_file_comments --cpu 10 -i $id -o $id);done
cat 6.Function_annotation/eggNOG/*.chunk_*.emapper.seed_orthologs > 6.Function_annotation/eggNOG/input_file.emapper.seed_orthologs
/media/ruizhi/software/eggnog-mapper/emapper.py --annotate_hits_table 6.Function_annotation/eggNOG/input_file.emapper.seed_orthologs --no_file_comments -o 6.Function_annotation/eggNOG/output_file --cpu 10
echo "query_name\tseed eggNOG ortholog\tseed ortholog evalue\tseed ortholog score\tPredicted taxonomic group\tPredicted protein name\tGene Ontology terms \tEC number\tKEGG_ko\tKEGG_Pathway\t KEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy \tBiGG Reaction\ttax_scope: eggNOG taxonomic level used for annotation\teggNOG OGs \tbestOG (deprecated, use smallest from eggnog OGs)\tCOG Functional Category\teggNOG free text description" > 6.Function_annotation/eggNOG/head.xls
cat 6.Function_annotation/eggNOG/head.xls 6.Function_annotation/eggNOG/output_file.emapper.annotations > 6.Function_annotation/eggNOG/output_file.emapper.annotations.xls
#GO
mkdir -p 6.Function_annotation/GO/
python2 /media/ruizhi/software/GO_annotation/GO.py 6.Function_annotation/eggNOG/output_file.emapper.annotations > 6.Function_annotation/GO/GO.annotation
cd 6.Function_annotation/GO/ 
python3 /media/ruizhi/software/GO_annotation/GO.ont.py -i GO.annotation -o GO.annotation.stast
Rscript /media/ruizhi/software/GO_annotation/barplot.r GO.secondary.annotation.csv
python before_GO_stast.py -i GO.annotation -o GO_forstat


cd ../../

#KEGG
mkdir -p 6.Function_annotation/KEGG/
split -l 200000 -a 3 -d 5.Non_redundunt_geneset/non_redundunt.nt.format.fa 6.Function_annotation/KEGG/non_redundunt.format.fa_chunk_
cd 6.Function_annotation/KEGG/
ls non_redundunt.format.fa_chunk_* | grep -v *.m8 | while read id;do(diamond blastp -p 10 -d /media/ruizhi/database/KEGG/kegg/genes/fasta/prokaryotes -e 1e-5 -k 1 --sensitive -q $id -o $(basename $id).kegg.m8);done
cat *.m8 > non_redundunt.nt.format.fa.kegg.m8
echo -e "Query_id\tSubject_id\tIdentity\tAlign_length\tMismatch\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE-value\tScore" > head.xls
cat head.xls non_redundunt.nt.format.fa.kegg.m8 > non_redundunt.nt.format.fa.kegg.xls
python /media/ruizhi/database/KEGG/kegg/genes/fasta/script/ko_id.py -k /media/ruizhi/database/KEGG/kegg/genes/fasta/script/KO_id_dict -p /media/ruizhi/database/KEGG/kegg/genes/fasta/script/path_id_dict -a /media/ruizhi/database/KEGG/kegg/genes/fasta/script/ko_path_id_dict -b 6.Function_annotation/KEGG/non_redundunt.nt.format.fa.kegg.m8 -o 6.Function_annotation/KEGG/non_redundunt.nt.kegg.annotation
cd ../../
cd ..
######################################################################################################
#ls *R1.fq.gz | while read id;do(cd-hit -i Results/6.Function_annotation/$(basename $id "_R1.fq.gz").faa -o Results/5.Non_redundunt_geneset/$(basename $id "_R1.fq.gz").non_redundunt.pr.fa);done

mkdir -p Results/7.Taxonomy/metaphlan
mv ~/miniconda3/bin/metaphlan_databases/mpa_previous ~/miniconda3/bin/metaphlan_databases/mpa_latest
#ls *.unmapped.R1.fq | while read id;do(cat $id $(basename $id '.unmapped.R1.fq').unmapped.R2.fq > $(basename $id '.unmapped.R1.fq').merge.fq);done
ls *.merge.fq | while read id;do(metaphlan2.py $id --input_type fastq --nproc 10 > Results/7.Taxonomy/metaphlan/$(basename $id '.merge.fq')_profile.txt);done

merge_metaphlan_tables.py Results/7.Taxonomy/metaphlan/*_profile.txt | sed 's/metaphlan2_//g' > Results/7.Taxonomy/metaphlan/merged_metaphlan2.txt

cat Results/7.Taxonomy/metaphlan/merged_metaphlan2.txt | grep s__ > Results/7.Taxonomy/metaphlan/1.species.merged_metaphlan2.txt
awk -F "__" '{print $NF}' Results/7.Taxonomy/metaphlan/1.species.merged_metaphlan2.txt > Results/7.Taxonomy/metaphlan/2.species.merged_metaphlan2.txt
head -n 2 Results/7.Taxonomy/metaphlan/merged_metaphlan2.txt > Results/7.Taxonomy/metaphlan/head.txt
sed -i '1d' Results/7.Taxonomy/metaphlan/head.txt
cat Results/7.Taxonomy/metaphlan/head.txt Results/7.Taxonomy/metaphlan/2.species.merged_metaphlan2.txt > Results/7.Taxonomy/metaphlan/species.merged_metaphlan2.txt
sed -i 's/_profile//g' Results/7.Taxonomy/metaphlan/species.merged_metaphlan2.txt
rm Results/7.Taxonomy/metaphlan/head.txt Results/7.Taxonomy/metaphlan/2.species.merged_metaphlan2.txt Results/7.Taxonomy/metaphlan/1.species.merged_metaphlan2.txt
hclust2.py -i  Results/7.Taxonomy/metaphlan/species.merged_metaphlan2.txt -o  Results/7.Taxonomy/metaphlan/abundance_heatmap_species.png --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
cd  Results/7.Taxonomy/metaphlan/
Rscript /media/ruizhi/software/metaphlan.abundunce.r species.merged_metaphlan2.txt


python2 /media/ruizhi/software/graphlan/export2graphlan/export2graphlan.py --skip_rows 1,2 -i merged_metaphlan2.txt --tree merged_metaphlan2.tree.txt --annotation merged_metaphlan2.annot.txt --most_abundant 40 --abundance_threshold 1 --least_biomarkers 10 --annotations 3,4,5,6 --external_annotations 7 --min_clade_size 1

# 绘图
conda activate python27
graphlan_annotate.py --annot merged_metaphlan2.annot.txt merged_metaphlan2.tree.txt merged_metaphlan2.abundance.xml
conda activate base
python2 /media/ruizhi/software/graphlan/graphlan.py --dpi 300 merged_metaphlan2.abundance.xml merged_abundance.png --external_legends
cd ../../

#######################################################################################################
#PCA  PCOA  normal_tree Alpha
mkdir -p Results/kraken
mkdir -p Results/bracken

ls *.unmapped.R1.fq | while read id;do(kraken --threads 20 --db /media/ruizhi/software/kraken_db/minikraken_20171019_8GB --preload --paired --fastq-input $id $(basename $id 'unmapped.R1.fq')unmapped.R2.fq --output Results/kraken/$(basename $id 'unmapped.R1.fq')txt;kraken-report --db /media/ruizhi/software/kraken_db/minikraken_20171019_8GB Results/kraken/$(basename $id 'unmapped.R1.fq')txt > Results/kraken/$(basename $id 'unmapped.R1.fq')kraken.report);done
#ls *.unmapped.R1.fq | while read id;do(kraken2 --db /media/ruizhi/software/kraken_db/minikraken_20171019_8GB  --threads 20  --report Results/kraken/$(basename $id 'unmapped.R1.fq')kraken.report --output Results/kraken/$(basename $id 'unmapped.R1.fq')txt  --paired $id $(basename $id 'unmapped.R1.fq')unmapped.R2.fq);done
ls *.unmapped.R1.fq | while read id;do(bracken -d /media/ruizhi/software/kraken_db/minikraken_20171019_8GB -i Results/kraken/$(basename $id 'unmapped.R1.fq')kraken.report -o Results/bracken/$(basename $id 'unmapped.R1.fq')bracken.report);done

cd Results/kraken
ls *.kraken_bracken.report | while read id;do(kreport2mpa.py -r $id -o ../../Results/bracken/$(basename $id 'kraken_bracken.report')txt);done
cd ../../
cd Results/bracken

########
ls *.txt | while read id;do(mv $id $(basename $id "txt")csv);done
python2 /media/ruizhi/software/Bracken/src/combine_bracken_outputs.py --files *.report -o total.bracken.file
sed -i 's/bracken.report.o_//g' total.bracken.file

sed 's/\t/,/g' total.bracken.file > test
Rscript /media/ruizhi/software/tran.r  test
sed -i 's/"//g' t_output
sed -i '2d' t_output
cp t_output species.abundunce.xls
rm test
awk 'NR%2==1{print}' t_output > t_output2
Rscript /media/ruizhi/software/tran.r t_output2
sed -i 's/"//g' t_output
sed -i 's/.num//g' t_output
rm t_output2


cat *.csv > total.test
cut -f 1 total.test > total.test2
cat total.test2 | grep s_ > total.test3
sort total.test3 | uniq > total.test4
python /media/ruizhi/software/anno.py -i total.test4 -o map_file -d t_output -m t.map_output
head -n 1 map_file > head.map_file
sed '1d' map_file | sort > sort.map_file
cat head.map_file sort.map_file > sort.head.map_file
head -n 1 t.map_output > head.t.map_output
sed '1d' t.map_output | sort > sort.t.map_output
cat head.t.map_output sort.t.map_output > sort.head.t.map_output

rm total.test* map_file head.map_file sort.map_file t.map_output head.t.map_output sort.t.map_output t_output
mv sort.head.map_file species.annotation.xls
mv total.bracken.file species.num.frac.xls
mv sort.head.t.map_output species.abundunce.xls
cd ../
mkdir 7.Taxonomy/
mv bracken/ 7.Taxonomy/
mv kraken/ 7.Taxonomy/
mv metaphlan/ 7.Taxonomy/


#generate sort.head.t.map_output and sort.head.map_file
#######################################################################################################


#build index

conda activate python37
cd Results/5.Non_redundunt_geneset/
ls non_redundunt.fnn.format.fa | while read id;do(salmon index -t $id -i transcript_index --type quasi -k 31 --perfectHash);done
cd ../../

#for file in *counts
#do
  # 提取样品名
#  name=${file%%.*}
  # 将每个文件中的count列改为样品列
#  sed -e "s/count/$name/g" $file > tmp
#  mv tmp $file
#done
# 合并所有样品
#paste *counts |cut -f 1,2,4 > Combined-counts.tab

#gene_abundunce
mkdir -p Results/8.Gene_abundunce/
ls *.unmapped.R1.fq | while read id;do(salmon quant -i Results/5.Non_redundunt_geneset/transcript_index --libType IU -1 $id -2 $(basename $id ".unmapped.R1.fq").unmapped.R2.fq -o Results/8.Gene_abundunce/$(basename $id '.unmapped.R1.fq'));done
ls *.merge.fq | while read id;do(name=$(basename $id '.merge.fq')_TPM;sed -e "s/TPM/$name/g" Results/8.Gene_abundunce/$(basename $id '.merge.fq')/quant.sf > Results/8.Gene_abundunce/$(basename $id '.merge.fq')/$(basename $id '.merge.fq').quant.sf);done
cut -f 1 Results/8.Gene_abundunce/*/*quant.sf > Results/8.Gene_abundunce/colname
ls *.merge.fq | while read id;do(cut -f 4 Results/8.Gene_abundunce/$(basename $id '.merge.fq')/$(basename $id '.merge.fq').quant.sf > Results/8.Gene_abundunce/$(basename $id '.merge.fq')/$(basename $id '.merge.fq').tmp);done

paste Results/8.Gene_abundunce/colname Results/8.Gene_abundunce/*/*.tmp > Results/8.Gene_abundunce/total.TPM
#need col.names


ls *.merge.fq | while read id;do(name=$(basename $id '.merge.fq')_NumReads;sed -e "s/NumReads/$name/g" Results/8.Gene_abundunce/$(basename $id '.merge.fq')/quant.sf > Results/8.Gene_abundunce/$(basename $id '.merge.fq')/$(basename $id '.merge.fq').NR.quant.sf);done

ls *.merge.fq | while read id;do(cut -f 5 Results/8.Gene_abundunce/$(basename $id '.merge.fq')/$(basename $id '.merge.fq').NR.quant.sf > Results/8.Gene_abundunce/$(basename $id '.merge.fq')/$(basename $id '.merge.fq').tmp);done

paste Results/8.Gene_abundunce/colname Results/8.Gene_abundunce/*/*.tmp > Results/8.Gene_abundunce/total.NumReads
#need col.names
rm  Results/8.Gene_abundunce/colname 
conda activate base

mkdir 9.Differential analysis
cd ../../Results/5.Non_redundunt_geneset/6.Function_annotation/
python /media/ruizhi/software/GO_annotation/before_GO_stast.py -i GO.annotation -o GO_forstat
Rscript /media/ruizhi/software/GO__annotation/GO_enrich.r GO_forstat

python /media/ruizhi/software/KEGG_annotation/KEGG.stast.py -i output_file.emapper.annotations -k kegg -l kegg.list -o kegg.ko.lis
#kegg.ko.lis for KEGG.enrich.r

Rscript /media/ruizhi/software/KEGG_annotation/KEGG.enrich.r kegg.ko.list /media/ruizhi/NGS/20210.7.16.metagenomes/Cleandata/Results/8.Gene_abundunce/test.diff_gene_deseq2.csv





