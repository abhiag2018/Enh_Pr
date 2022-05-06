#!/bin/bash
#SBATCH --mem=100G
#SBATCH --qos batch
#SBATCH --time 1-00:00:00

species=$1 ## "hg" # hg / mm
max_distance=$2 ## 1e6
enhancer_window=$3 ## 2000
promoter_window=$4 ## 1000
enh_inp=$5 # "enhancer.bed"
pr_inp=$6 # "promoter.bed"

if [[ ${species:0:2} == "hg" ]]
then
	chrom_list=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
	echo -e "processing all human chromosomes \n"
else if [[ ${species:0:2} == "mm" ]]
	then
		chrom_list=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
		echo -e "processing all mouse chromosomes \n"
	else
		echo -e "unknown species $species \n"
		exit 1
	fi
fi


awk 'BEGIN{OFS="\t"}{print $1,(($2+$3))/2,$4}' $enh_inp > enh.bed
awk 'BEGIN{OFS="\t"}{print $1,(($2+$3))/2,$4}' $pr_inp > pr.bed
num=0
for chrom in "${chrom_list[@]}"
do
	echo "processing chromosome $chrom.."
	awk -v chrom="$chrom" -v OFS="\t" '{ if ( $1==chrom ) print $1,$2,$3}' pr.bed > pr_$chrom.bed 
	awk -v chrom="$chrom" -v OFS="\t" '{ if ( $1==chrom ) print $1,$2,$3}' enh.bed > enh_$chrom.bed 

	python -c "import pandas as pd; \
	df_pair1=pd.read_csv('enh_$chrom.bed',sep='\t',header=None).merge(pd.read_csv('pr_$chrom.bed',sep='\t',header=None),how='cross'); \
	df_pair=df_pair1[abs(df_pair1['1_x']-df_pair1['1_y'])<$max_distance].copy(deep=True); \
	df_pair['1_x'] = df_pair['1_x']-$enhancer_window//2; \
	df_pair['1_y'] = df_pair['1_y']-$promoter_window//2; \
	df_pair=df_pair.rename(columns={'0_x':'enhancer_chrom','1_x':'enhancer_start','2_x':'enhancer_name','0_y':'promoter_chrom','1_y':'promoter_start','2_y':'promoter_name'}); \
	df_pair['enhancer_end']=df_pair['enhancer_start']+$enhancer_window; \
	df_pair['promoter_end']=df_pair['promoter_start']+$promoter_window; \
	df_pair['label']=1; \
	col_list = ['enhancer_chrom','enhancer_start','enhancer_end','enhancer_name','promoter_chrom','promoter_start','promoter_end','promoter_name','label']; \
	df_pair[col_list].to_csv('$chrom.csv',index=False)"
	rm pr_$chrom.bed enh_$chrom.bed
	num=$(($num+`wc -l $chrom.csv | awk '{print $1}'`-1))
done

cp chr1.csv pairs_bruteForce.$species.csv
for chrom in "${chrom_list[@]:1}"
do
	tail -n +2 $chrom.csv >> pairs_bruteForce.$species.csv
done

# rm chr*.csv enh.bed pr.bed
