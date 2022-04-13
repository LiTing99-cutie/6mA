################################################
#File Name: /home/user/data2/lit/project/6mA/bin/binned_level.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 12 Apr 2022 04:13:04 PM CST
################################################

#!/bin/sh 

set -eou pipefail

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
  echo -e "     -fs|--fastaSize [Mandatory] \t\t\t\tgenome sorted chromosome size"
  echo -e "     -fa|--fasta [Mandatory] \t\t\t\tgenome fasta"
  echo -e "     -ld|--logDir [mandatory] \t\t\t\tlog directory"
  echo -e "     -bd|--binned [mandatory] \t\t\t\binned bed directory"
  echo -e "     -h|--help [optional] \t\t\t\tprint this help page"
}

## Default argument

while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--sample)             sample=$2;shift;;
        -mA|--mA)           mA_list=$2;shift;;
        -o|--ouput)           output_dir=$2;shift;;
        -fs|--fastaSize)           genome_sorted_chrSize=$2;shift;;
        -fa|--fasta)            genome_fasta=$2;shift;;
        -ld|--logDir)                logDir=$2;shift;;
        -bd|--binned)               binned_dir=$2;shift;;    
        -h)                                     usage;exit 1;;
        -|--)                                     shift; break;;
        *)                                      usage; echo -e "\n[ERR] $(date) Unkonwn option: $1"; exit 1;;
   esac
    shift
done

## Check mandatory argument
[[ -f $genome_sorted_chrSize  ]] || { usage;echo -e "\n[ERR] $(date) -fs not assigned";exit; }
[[ -f $genome_fasta  ]] ||  { usage;echo -e "\n[ERR] $(date) -fa not assigned";exit; }
[[ -n $sample  ]] || { usage;echo -e "\n[ERR] $(date) -s not assigned";exit; }
[[ -n $mA_list  ]] || { usage;echo -e "\n[ERR] $(date) -mA not assigned";exit; }
[[ -n $output_dir  ]] || { usage;echo -e "\n[ERR] $(date) -o not assigned";exit; }
[[ -n $logDir  ]] || { usage;echo -e "\n[ERR] $(date) -ld not assigned";exit; }
[[ -n $binned_dir  ]] || { usage;echo -e "\n[ERR] $(date) -bd not assigned";exit; }
[[ -d $output_dir ]] || mkdir -p $output_dir
[[ -d $logDir ]] || mkdir -p $logDir
[[ -f $logDir/run.log ]] && rm -rf $logDir/run.log

echo -e "sample\t$sample"
echo -e "mA_list\t$mA_list"
echo -e "output_dir\t$output_dir"
echo -e "genome_sorted_chrSize\t$genome_sorted_chrSize"
echo -e "genome_fasta\t$genome_fasta"
echo -e "logDir\t$logDir"
echo -e "binned_dir\t$binned_dir"

SAMPLE=$( echo $sample | awk -F"|" -v OFS=" " '{for(i=1;i<=NF;i++) print $i}' )
MA_LIST=$( echo $mA_list | awk -F"|" -v OFS=" " '{for(i=1;i<=NF;i++) print $i}' )

    ## 2.2 6mA level ###

        # mA level
        time (
        for sample in ${SAMPLE};do
            for mA_list in ${MA_LIST};do
                for j in nucleosome linker;do
                    for i in `seq 1 10`;do
                        BASENAME=`basename $mA_list`
                        list_tmp=${mA_list%%.*} && list=${list_tmp##*_} 
                        echo "*** process ${sample}.$j.${list}.bin_$i"
                        less $binned_dir/${sample}.$j.binned.bed | \
                                grep $i$ | sort -k1,1 -k2,2n > \
                                    $output_dir/${sample}.$j.bin_$i.bed
                                
                        less $output_dir/${sample}.$j.bin_$i.bed | awk -v OFS='\t' '{if($1=="MtDNA") print "chrM",$2,$3,$4 ; else print "chr"$1,$2,$3,$4}' | \
                        sort -k1,1 -k2,2n > $output_dir/${sample}.$j.bin_$i.chr.bed

                        bedtools intersect -wa -wb -a $output_dir/${sample}.$j.bin_$i.chr.bed -b ${mA_list} -sorted > $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed &
                    done
                    wait
                done 
            done
        done ) >> $logDir/run.log 2>&1

        time (
        for sample in ${SAMPLE};do
            for j in nucleosome linker;do
                for i in `seq 1 10`;do
                    echo "*** process ${sample}.$j.bin_$i"
                    bedtools nuc -fi $genome_fasta -bed $output_dir/${sample}.$j.bin_$i.bed > \
                    $output_dir/${sample}.$j.bin_$i.txt &
                done
                wait
            done 
        done ) >> $logDir/run.log 2>&1

        [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt

        for sample in ${SAMPLE};do
            for mA_list in ${MA_LIST};do
                for j in nucleosome linker;do
                    for i in `seq 1 10`;do
                        BASENAME=`basename $mA_list`
                        list_tmp=${mA_list%%.*} && list=${list_tmp##*_}
                        mA=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | wc -l`
                        A=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' `
                        C=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' `
                        G=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' `
                        T=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' `
                        base_number=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $13} END {print sum}' `
                        AT=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' `
                        GC=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' `
                        FreqSum=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $8} END {print sum}'`
                        echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${C}\t${G}\t${T}\t${base_number}\t${AT}\t${GC}\t${FreqSum}" >> $output_dir/binned.level.txt
                    done
                done
            done
        done 
    ### EO 2.2 6mA level ###

### EO 2. binned linker or nucleosome enrichment ###


### merge different samples 

    mkdir -p $output_dir/merge   

    for sample in ${SAMPLE};do
        for mA_list in ${MA_LIST};do
            for j in nucleosome;do
                BASENAME=`basename $mA_list`
                list_tmp=${mA_list%%.*} && list=${list_tmp##*_}
                cat $output_dir/${sample}.$j.${list}.bin_5.intersect.bed $output_dir/${sample}.$j.${list}.bin_6.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_1.intersect.bed
                cat $output_dir/${sample}.$j.${list}.bin_4.intersect.bed $output_dir/${sample}.$j.${list}.bin_7.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_2.intersect.bed
                cat $output_dir/${sample}.$j.${list}.bin_3.intersect.bed $output_dir/${sample}.$j.${list}.bin_8.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_3.intersect.bed
                cat $output_dir/${sample}.$j.${list}.bin_2.intersect.bed $output_dir/${sample}.$j.${list}.bin_9.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_4.intersect.bed
                cat $output_dir/${sample}.$j.${list}.bin_1.intersect.bed $output_dir/${sample}.$j.${list}.bin_10.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_5.intersect.bed
                cat $output_dir/${sample}.$j.bin_5.txt $output_dir/${sample}.$j.bin_6.txt > $output_dir/merge/${sample}.$j.bin_1.txt
                cat $output_dir/${sample}.$j.bin_4.txt $output_dir/${sample}.$j.bin_7.txt > $output_dir/merge/${sample}.$j.bin_2.txt
                cat $output_dir/${sample}.$j.bin_3.txt $output_dir/${sample}.$j.bin_8.txt > $output_dir/merge/${sample}.$j.bin_3.txt
                cat $output_dir/${sample}.$j.bin_2.txt $output_dir/${sample}.$j.bin_9.txt > $output_dir/merge/${sample}.$j.bin_4.txt
                cat $output_dir/${sample}.$j.bin_1.txt $output_dir/${sample}.$j.bin_10.txt > $output_dir/merge/${sample}.$j.bin_5.txt
            done    
        done
    done

    for sample in ${SAMPLE};do
        for mA_list in ${MA_LIST};do
            for j in linker;do
                BASENAME=`basename $mA_list`
                list_tmp=${mA_list%%.*} && list=${list_tmp##*_}
                cat $output_dir/${sample}.$j.${list}.bin_5.intersect.bed $output_dir/${sample}.$j.${list}.bin_6.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_5.intersect.bed
                cat $output_dir/${sample}.$j.${list}.bin_4.intersect.bed $output_dir/${sample}.$j.${list}.bin_7.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_4.intersect.bed
                cat $output_dir/${sample}.$j.${list}.bin_3.intersect.bed $output_dir/${sample}.$j.${list}.bin_8.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_3.intersect.bed
                cat $output_dir/${sample}.$j.${list}.bin_2.intersect.bed $output_dir/${sample}.$j.${list}.bin_9.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_2.intersect.bed
                cat $output_dir/${sample}.$j.${list}.bin_1.intersect.bed $output_dir/${sample}.$j.${list}.bin_10.intersect.bed > \
                $output_dir/merge/${sample}.$j.${list}.bin_1.intersect.bed
                cat $output_dir/${sample}.$j.bin_5.txt $output_dir/${sample}.$j.bin_6.txt > $output_dir/merge/${sample}.$j.bin_5.txt
                cat $output_dir/${sample}.$j.bin_4.txt $output_dir/${sample}.$j.bin_7.txt > $output_dir/merge/${sample}.$j.bin_4.txt
                cat $output_dir/${sample}.$j.bin_3.txt $output_dir/${sample}.$j.bin_8.txt > $output_dir/merge/${sample}.$j.bin_3.txt
                cat $output_dir/${sample}.$j.bin_2.txt $output_dir/${sample}.$j.bin_9.txt > $output_dir/merge/${sample}.$j.bin_2.txt
                cat $output_dir/${sample}.$j.bin_1.txt $output_dir/${sample}.$j.bin_10.txt > $output_dir/merge/${sample}.$j.bin_1.txt
            done    
        done
    done

    merge_output_dir=$output_dir/merge
    [ -f $merge_output_dir/binned.level.txt ] && rm -rf $merge_output_dir/binned.level.txt
    for sample in ${SAMPLE};do
        for mA_list in ${MA_LIST};do
            for j in nucleosome linker;do
                for i in `seq 1 5`;do
                    BASENAME=`basename $mA_list`
                    list_tmp=${mA_list%%.*} && list=${list_tmp##*_}
                    mA=`less $merge_output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | wc -l`
                    A=` sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' `
                    C=` sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' `
                    G=` sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' `
                    T=` sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' `
                    base_number=` sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $13} END {print sum}' `
                    AT=` sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' `
                    GC=` sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' `
                    FreqSum=`less $merge_output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $8} END {print sum}'`
                    echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${C}\t${G}\t${T}\t${base_number}\t${AT}\t${GC}\t${FreqSum}" >> $merge_output_dir/binned.level.txt
                done
            done 
        done
    done
### EO　merge different samples 
