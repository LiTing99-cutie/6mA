################################################
#File Name: /home/user/data2/lit/project/6mA/bin/binned.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 12 Apr 2022 05:47:26 PM CST
################################################

#!/bin/sh 

set -eou pipefail

usage(){
  echo "Usage: bash $(basename $0) -s sample -mA mA_list -o output_dir -c cutoff -fs genome_sorted_chrSize -fa genome_fasta -ld logDir [-h]"
  echo "Author: LiT"
  echo "Description: This script binned nucleosome and linker in whole genome and generate binned.bed in correspond directory."
  echo "Date: 2022-03-23"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     -s|--sample [Mandatory] \t\t\t\tsample"
  echo -e "     -o|--output_dir [Mandatory] \t\t\t\toutput directory"
  echo -e "     -fs|--fastaSize [Mandatory] \t\t\t\tgenome sorted chromosome size"
  echo -e "     -ld|--logDir [mandatory] \t\t\t\tlog directory"
  echo -e "     -h|--help [optional] \t\t\t\tprint this help page"
}

## Default argument

while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--sample)             sample=$2;shift;;
        -o|--ouput)           output_dir=$2;shift;;
        -fs|--fastaSize)           genome_sorted_chrSize=$2;shift;;
        -ld|--logDir)                logDir=$2;shift;;
        -h)                                     usage;exit 1;;
        -|--)                                     shift; break;;
        *)                                      usage; echo -e "\n[ERR] $(date) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

## Check mandatory argument
[[ -f $genome_sorted_chrSize  ]] || { usage;echo -e "\n[ERR] $(date) -fs not assigned";exit; }
[[ -n $sample  ]] || { usage;echo -e "\n[ERR] $(date) -s not assigned";exit; }
[[ -n $output_dir  ]] || { usage;echo -e "\n[ERR] $(date) -o not assigned";exit; }
[[ -n $logDir  ]] || { usage;echo -e "\n[ERR] $(date) -ld not assigned";exit; }
[[ -d $output_dir ]] || mkdir -p $output_dir
[[ -d $logDir ]] || mkdir -p $logDir
[[ -f $logDir/run.log ]] && rm -rf $logDir/run.log

echo -e "sample\t$sample"
echo -e "output_dir\t$output_dir"
echo -e "genome_sorted_chrSize\t$genome_sorted_chrSize"
echo -e "logDir\t$logDir"

SAMPLE=$( echo $sample | awk -F"|" -v OFS=" " '{for(i=1;i<=NF;i++) print $i}' )

### 1 get bed ###

### 1.1 nucleosome bed ###
  time (
  for sample in ${SAMPLE};do
      echo "*** generating ${sample}.nucleosome.bed"
      # modify here
      sed '1d' danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls | \
      awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $output_dir/${sample}.nucleosome.bed 
      nucleosome_number=`less $output_dir/${sample}.nucleosome.bed | wc -l`
      echo "${sample} nucleosome number: ${nucleosome_number}"
  done ) >> $logDir/binned.log 2>&1
### EO 1.1 nucleosome bed ###

### 1.2 linker bed ###
  time (
  for sample in ${SAMPLE};do
        echo "*** generating ${sample}.linker.bed"
        # do not remove the first linker and the last linker
        bedtools complement  -i $output_dir/${sample}.nucleosome.bed -g $genome_sorted_chrSize > \
        $output_dir/${sample}.linker.bed

        linker_number=`less $output_dir/${sample}.linker.bed | wc -l`
        echo "${sample} linker_number: ${linker_number}"
  done ) >> $logDir/binned.log 2>&1
### EO 1.2 linker bed ###

### EO 1. get bed ###

### 2. binned linker or nucleosome enrichment ###

    ### 2.1 bin ### 
        time (
        for sample in ${SAMPLE};do
            echo "bin ${sample}"
            bedtools makewindows -b $output_dir/${sample}.nucleosome.bed  -n 10 -i winnum > \
            $output_dir/${sample}.nucleosome.binned.bed  2> $logDir/${sample}.nucleosome.bin.log &
            bedtools makewindows -b $output_dir/${sample}.linker.bed  -n 10 -i winnum > \
            $output_dir/${sample}.linker.binned.bed 2> $logDir/${sample}.linker.bin.log &
        done 
        wait ) >> $logDir/binned.log 2>&1
    ### EO 2.1 bin ###