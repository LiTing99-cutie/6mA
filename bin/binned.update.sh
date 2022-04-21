################################################
#File Name: /home/user/data2/lit/project/6mA/bin/binned.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 12 Apr 2022 05:47:26 PM CST
################################################

#!/bin/sh 

set -eou pipefail

sample=$1 # Mnase sample
output_dir=$2 #output directory
genome_sorted_chrSize=$3 # genome sorted chromosome size, sort -k1,1 -k2,2n 
nucleosome_bed=$4 # no header, chr start end ...
centerORnot=$5 # binned directly or binned conservely
bin_number=$6 # number of bins

echo "sample: $sample"
echo "output_dir: $output_dir"
echo "genome_sorted_chrSize: $genome_sorted_chrSize"
echo "nucleosome_bed: $nucleosome_bed"
echo "centerORnot: $centerORnot"

[ -d $output_dir ] || mkdir -p $output_dir

[ -f ${genome_sorted_chrSize}.bed ] ||  { echo -e "\n[ERR] $(date) ${genome_sorted_chrSize}.bed not exist";exit; }

if [ ! centerORnot ] ;then
  # 1 nucleosome bed
    time (
    echo "*** not center: generating ${sample}.nucleosome.bed"
    less $nucleosome_bed | \
    awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $output_dir/${sample}.nucleosome.bed 
    nucleosome_number=`less $output_dir/${sample}.nucleosome.bed | wc -l`
    echo "${sample} nucleosome number: ${nucleosome_number}" 
    ) 
fi

if [ centerORnot ] ;then
  # 1 nucleosome bed
    time (
    echo "*** center 73 bp: generating ${sample}.nucleosome.bed"
    less $nucleosome_bed | \
    awk -v OFS="\t" '{print $1,(($2+$3+10)/2-73),(($2+$3+10)/2+73)}' | sort -k1,1 -k2,2n | grep -v '-' | \
    bedtools intersect -a - -b ${genome_sorted_chrSize}.bed  > $output_dir/${sample}.nucleosome.bed 
    nucleosome_number=`less $output_dir/${sample}.nucleosome.bed | wc -l`
    echo "${sample} nucleosome number: ${nucleosome_number}" 
    ) 
fi


# 2 linker bed 
  time (
  echo "*** generating ${sample}.linker.bed"
  bedtools complement  -i $output_dir/${sample}.nucleosome.bed -g $genome_sorted_chrSize > \
  $output_dir/${sample}.linker.bed

  linker_number=`less $output_dir/${sample}.linker.bed | wc -l`
  echo "${sample} linker_number: ${linker_number}" ) 

# 3.bin
  time (
    echo "bin ${sample}"
    bedtools makewindows -b $output_dir/${sample}.nucleosome.bed  -n $bin_number -i winnum > \
    $output_dir/${sample}.nucleosome.binned.bed  2>/dev/null & 
    bedtools makewindows -b $output_dir/${sample}.linker.bed  -n $bin_number -i winnum > \
    $output_dir/${sample}.linker.binned.bed 2>/dev/null &
  wait ) 

# 4. binned bed
for j in nucleosome linker;do
  for i in `seq 1 $bin_number`;do
    echo "*** process ${sample}.$j.bin_$i"
    less $output_dir/${sample}.$j.binned.bed | \
    awk -v OFS='\t' '{if ($4=="'$i'") print $0}' > \
                $output_dir/${sample}.$j.bin_$i.bed
  done
done