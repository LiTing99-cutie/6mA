################################################
#File Name: /home/user/data2/lit/project/6mA/bin/filter_fuzziness_score_bin_level.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 23 Mar 2022 11:13:45 AM CST
################################################

#!/bin/sh 

set -e

usage(){
  echo "Usage: bash $(basename $0) -s sample -mA mA_list -o output_dir -c cutoff -fs genome_sorted_chrSize -fa genome_fasta -ld logDir [-h]"
  echo "Author: LiT"
  echo "Description: This script calculate mA level and GC/AT level in binned nucleosome and linker."
  echo "Date: 2022-03-23"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     -s|--sample [Mandatory] \t\t\t\tsample"
  echo -e "     -mA|--mA [Mandatory] \t\t\t\tmA list"
  echo -e "     -o|--output_dir [Mandatory] \t\t\t\toutput directory"
  echo -e "     -c|--cutoff [Mandatory] \t\t\t\tfuzziness score cutoff"
  echo -e "     -fs|--fastaSize [Mandatory] \t\t\t\tgenome sorted chromosome size"
  echo -e "     -fa|--fasta [Mandatory] \t\t\t\tgenome fasta"
  echo -e "     -ld|--logDir [mandatory] \t\t\t\tlog directory"
  echo -e "     -h|--help [optional] \t\t\t\tprint this help page"
}

## Default argument

while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--sample)             sample=$2;shift;;
        -mA|--mA)           mA_list=$2;shift;;
        -o|--ouput)           output_dir=$2;shift;;
        -c|--cutoff)           cutoff=$2;shift;;
        -fs|--fastaSize)           genome_sorted_chrSize=$2;shift;;
        -fa|--fasta)            genome_fasta=$2;shift;;
        -ld|--logDir)                logDir=$2;shift;;
        -h)                                     usage;exit 1;;
        -|--)                                     shift; break;;
        *)                                      usage; echo -e "\n[ERR] $(date) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

## Check mandatory argument
echo $genome_sorted_chrSize
[[ -f $genome_sorted_chrSize  ]] || { usage;echo -e "\n[ERR] $(date) -fs not assigned";exit; }
[[ -f $genome_fasta  ]] ||  { usage;echo -e "\n[ERR] $(date) -fa not assigned";exit; }
[[ -n $sample  ]] || { usage;echo -e "\n[ERR] $(date) -s not assigned";exit; }
[[ -n $mA_list  ]] || { usage;echo -e "\n[ERR] $(date) -mA not assigned";exit; }
[[ -n $cutoff  ]] || { usage;echo -e "\n[ERR] $(date) -c not assigned";exit; }
[[ -n $output_dir  ]] || { usage;echo -e "\n[ERR] $(date) -o not assigned";exit; }
[[ -n $logDir  ]] || { usage;echo -e "\n[ERR] $(date) -ld not assigned";exit; }
[[ -d $output_dir ]] || mkdir -p $output_dir
[[ -d $logDir ]] || mkdir -p $logDir
[[ -f $logDir/run.log ]] && rm -rf $logDir/run.log

echo -e "sample\t$sample"
echo -e "mA_list\t$mA_list"
echo -e "output_dir\t$output_dir"
echo -e "cutoff\t$cutoff"
echo -e "genome_sorted_chrSize\t$genome_sorted_chrSize"
echo -e "genome_fasta\t$genome_fasta"
echo -e "logDir\t$logDir"

sample=$( echo $sample | awk -F"|" -v OFS=" " '{for(i=1;i<=NF;i++) print $i}' )
mA_list=$( echo $mA_list | awk -F"|" -v OFS=" " '{for(i=1;i<=NF;i++) print $i}' )

### 1 get bed ###

### 1.1 nucleosome bed ###
  time (
  for sample in ${sample};do
      echo "*** generating ${sample}.nucleosome.bed"
      # modify here
      sed '1d' danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls | \
      awk -v OFS="\t" '{print $1,$2,$3,NR}' | sort -k1,1 -k2,2n > $output_dir/${sample}.nucleosome.bed 
      # modify here
      sed '1d' danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls | \
      awk -v OFS='\t' '{print $0,NR}' | sort -k1,1 -k2,2n | \
      # filter based on fuzziness score
      awk -v OFS="\t" '{if($6<'$cutoff') print $1,$2,$3,$7}' > $output_dir/${sample}.nucleosome.filter.bed
      nucleosome_number=`less $output_dir/${sample}.nucleosome.filter.bed | wc -l`
      echo "${sample} nucleosome number after filtering based on fuzziness score(cutoff:$cutoff): ${nucleosome_number}"
  done ) >> $logDir/run.log 2>&1
### EO 1.1 nucleosome bed ###
set +e
### 1.2 linker bed ###
  time (
  for sample in ${sample};do
      echo "*** generating ${sample}.linker.bed"

      bedtools complement  -i $output_dir/${sample}.nucleosome.bed -g ${genome_sorted_chrSize} | awk -v OFS="\t" '{print $1,$2,$3,NR}' > \
      $output_dir/${sample}.linker.bed

      cut -f 4  $output_dir/${sample}.nucleosome.filter.bed | awk '{print $1+1}' > $output_dir/${sample}.nucleosome_id.add_1.txt
      cut -f 4  $output_dir/${sample}.nucleosome.filter.bed > $output_dir/${sample}.nucleosome_id.txt
      cat $output_dir/${sample}.nucleosome_id.add_1.txt $output_dir/${sample}.nucleosome_id.txt | sort -k1,1n | uniq > $output_dir/${sample}.linker_id.txt

      join -1 1 -2 4 $output_dir/${sample}.linker_id.txt  $output_dir/${sample}.linker.bed > $output_dir/${sample}.linker.filter.bed 

      linker_number=`less $output_dir/${sample}.linker.filter.bed | wc -l`
      echo "${sample} linker number after filtering based on fuzziness score(cutoff:$cutoff): ${linker_number}"
  done ) >> $logDir/run.log 2>&1
### EO 1.2 linker bed ###
set -e

### EO 1. get bed ###

### 2. binned linker or nucleosome enrichment ###

### 2.1 bin ###
    time (
    for sample in ${sample};do
    echo "bin ${sample}"
        bedtools makewindows -b $output_dir/${sample}.nucleosome.bed  -n 10 -i winnum > \
        $output_dir/${sample}.nucleosome.binned.bed  2> $logDir/${sample}.nucleosome.bin.log &
        bedtools makewindows -b $output_dir/${sample}.linker.bed  -n 10 -i winnum > \
        $output_dir/${sample}.linker.binned.bed 2> $logDir/${sample}.linker.bin.log &
    done 
    wait ) >> $logDir/run.log 2>&1
### EO 2.1 bin ###

### 2.2 6mA level ###
    
    time ( for sample in ${sample};do
        for j in nucleosome linker;do
            for i in `seq 1 10`;do
                echo "*** process ${sample}.$j.bin_$i"
                less $output_dir/${sample}.$j.binned.bed | \
                grep $i$ > \
                $output_dir/${sample}.$j.bin_$i.bed
                    
                less $output_dir/${sample}.$j.bin_$i.bed | awk -v OFS='\t' '{if($1=="MtDNA") print "chrM",$2,$3,$4 ; else print "chr"$1,$2,$3,$4}' > \
                $output_dir/${sample}.$j.bin_$i.chr.bed

                # modify here
                bedtools intersect -wo -a danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.wig.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/${sample}.$j.intersect.occu.bin_$i.bed &

                bedtools nuc -fi $genome_fasta -bed $output_dir/${sample}.$j.bin_$i.bed > \
                $output_dir/${sample}.$j.bin_$i.txt &
            done
            wait
        done
    done ) >> $logDir/run.log 2>&1 )

    # mA level
    time (
    for sample in ${sample};do
        for mA_list in ${mA_list};do
            for j in nucleosome linker;do
                for i in `seq 1 10`;do
                    list_tmp=${mA_list%.*} && list=${list_tmp##*_} 
                    echo "*** process ${sample}.$j.${list}.bin_$i"
                    bedtools intersect -wa -wb -a $output_dir/${sample}.$j.bin_$i.chr.bed -b ${mA_list} > $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed &
                    done
                    wait
            done 
        done
    done ) >> $logDir/run.log 2>&1

    [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt

    time (
    for sample in ${sample};do
        for mA_list in ${mA_list};do
            for j in nucleosome linker;do
                for i in `seq 1 10`;do
                    list_tmp=${mA_list%.*} && list=${list_tmp##*_}
                    mA=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | wc -l`
                    A=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' `
                    C=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' `
                    G=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' `
                    T=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' `
                    base_number=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $13} END {print sum}' `
                    AT=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' `
                    GC=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' `
                    FreqSum=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $8} END {print sum}'`
                    OccuSum=`awk '{sum += $5*$10} END {print sum}' $output_dir/${sample}.$j.intersect.occu.bin_$i.bed`
                    echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${C}\t${G}\t${T}\t${base_number}\t${AT}\t${GC}\t${FreqSum}\t${OccuSum}" >> $output_dir/binned.level.txt
                done
            done
        done
    done ) >> $logDir/run.log 2>&1
    ### EO 2.2 6mA level ###

    ### EO 2. binned linker or nucleosome enrichment ###