#!/bin/sh 

################################################
#File Name: /home/user/data2/lit/project/6mA/bin/binned_level_Pacbio.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 15 Apr 2022 09:07:21 PM CST
################################################

set -eou pipefail

usage(){
  echo "Usage: bash $(basename $0) -s sample -mA mA_list -o output_dir -c cutoff -fs genome_sorted_chrSize -fa genome_fasta -ld logDir -bd binned_bed_dir -n bin_number [-h]"
  echo "Author: LiT"
  echo "Description: This script calculate mA level and GC/AT level in binned nucleosome and linker, and modification are called by Pacbio SMRTlink pipeline."
  echo "Date: 2022-04-15"
  echo "------------------------------------------------"
  echo "OPTIONS"
  echo -e "     -s|--sample [Mandatory] \t\t\t\tsample"
  echo -e "     -mA|--mA [Mandatory] \t\t\t\tmA list: same sequence name as nucleosome bed, like both chr or no chr; chr start end frequency or ipdratio; file name _*.bed, and * will bed used as list name "
  echo -e "     -o|--output_dir [Mandatory] \t\t\t\toutput directory"
  echo -e "     -fs|--fastaSize [Mandatory] \t\t\t\tgenome sorted chromosome size"
  echo -e "     -fa|--fasta [Mandatory] \t\t\t\tgenome fasta"
  echo -e "     -ld|--logDir [mandatory] \t\t\t\tlog directory"
  echo -e "     -bd|--binned [mandatory] \t\t\t\tbinned bed directory, with bin_1.bed, bin_2.bed..., next hierarchy should be samples correspond to -s"
  echo -e "     -n|--bin_number [mandatory] \t\t\t\tbin number"
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
        -n|--n)               bin_number=$2;shift;; 
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
[[ -n $bin_number  ]] || { usage;echo -e "\n[ERR] $(date) -n not assigned";exit; }
[[ -d $output_dir ]] || mkdir -p $output_dir

echo -e "sample\t$sample"
echo -e "mA_list\t$mA_list"
echo -e "output_dir\t$output_dir"
echo -e "genome_sorted_chrSize\t$genome_sorted_chrSize"
echo -e "genome_fasta\t$genome_fasta"
echo -e "logDir\t$logDir"
echo -e "binned_dir\t$binned_dir"
echo -e "bin_number\t$bin_number"


SAMPLE=$( echo $sample | awk -F"|" -v OFS=" " '{for(i=1;i<=NF;i++) print $i}' )
MA_LIST=$( echo $mA_list | awk -F"|" -v OFS=" " '{for(i=1;i<=NF;i++) print $i}' )
half_bin_number=$(($bin_number/2))


    ## 2.2 6mA level ###

        # mA level
        time (
        for sample in ${SAMPLE};do
            for mA_list in ${MA_LIST};do
                for j in nucleosome linker;do
                    for i in $(seq 1 $bin_number);do
                        BASENAME=$(basename $mA_list)
                        list_tmp=${BASENAME%%.*} && list=${list_tmp##*_}
                        echo "*** process ${sample}.$j.${list}.bin_$i"

                        bedtools intersect -wa -wb -a $binned_dir/${sample}/${sample}.$j.bin_$i.bed -b ${mA_list} -sorted > $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed &
                    done
                    wait
                done 
            done
        done ) 


        time (
        for sample in ${SAMPLE};do
            for j in nucleosome linker;do
                for i in $(seq 1 $bin_number);do
                    echo "*** process ${sample}.$j.bin_$i"
                    bedtools nuc -fi $genome_fasta -bed $binned_dir/${sample}/${sample}.$j.bin_$i.bed > \
                    $output_dir/${sample}.$j.bin_$i.txt &
                done
                wait
            done 
        done ) 

        [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt

        for sample in ${SAMPLE};do
            for mA_list in ${MA_LIST};do
                for j in nucleosome linker;do
                    for i in $(seq 1 $bin_number);do
                        BASENAME=$(basename $mA_list)
                        list_tmp=${BASENAME%%.*} && list=${list_tmp##*_}
                        mA=$(less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | wc -l)
                        A=$( sed '1d' "$output_dir"/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' )
                        C=$( sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' )
                        G=$( sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' )
                        T=$( sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' )
                        base_number=$( sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $13} END {print sum}' )
                        AT=$( sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' )
                        GC=$( sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' )
                        FreqSum=$(less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $8} END {print sum}')
                        echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${C}\t${G}\t${T}\t${base_number}\t${AT}\t${GC}\t${FreqSum}" >> $output_dir/binned.level.txt
                    done
                done
            done
        done 
    ### EO 2.2 6mA level ###

### EO 2. binned linker or nucleosome enrichment ###


### merge different samples 

    mkdir -p $output_dir/merge   
<<'!'
    for sample in ${SAMPLE};do
        for mA_list in ${MA_LIST};do
            for j in nucleosome;do
                BASENAME=`basename $mA_list`
                list_tmp=${BASENAME%%.*} && list=${list_tmp##*_}
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
!

type=nucleosome
for sample in ${SAMPLE};do
    for mA_list in ${MA_LIST};do
        BASENAME=$(basename $mA_list)
        list_tmp=${BASENAME%%.*} && list=${list_tmp##*_}
        i=0 && j=1 && k=1
        for loop in $(seq 1 $half_bin_number);do
            echo "*** loop $loop"
            x=$(($half_bin_number-$i)) && y=$(($half_bin_number+$j)) && z=$(($k))
            cat $output_dir/${sample}.$type.${list}.bin_$x.intersect.bed $output_dir/${sample}.$type.${list}.bin_$y.intersect.bed > \
            $output_dir/merge/${sample}.$type.${list}.bin_$z.intersect.bed &
            cat $output_dir/${sample}.$type.bin_$x.txt $output_dir/${sample}.$type.bin_$y.txt > $output_dir/merge/${sample}.$type.bin_$z.txt &
            i=$(($i+1)) && j=$(($j+1)) && k=$(($k+1))
        done
    done
done
wait


<<'!'
    for sample in ${SAMPLE};do
        for mA_list in ${MA_LIST};do
            for j in linker;do
                BASENAME=`basename $mA_list`
                list_tmp=${BASENAME%%.*} && list=${list_tmp##*_}
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
!


type=linker
for sample in ${SAMPLE};do
    for mA_list in ${MA_LIST};do
        BASENAME=$(basename $mA_list)
        list_tmp=${BASENAME%%.*} && list=${list_tmp##*_}
        i=0 && j=1 && k=$half_bin_number
        for loop in $(seq 1 $half_bin_number);do
            echo "*** loop $loop"
            x=$(($half_bin_number-$i)) && y=$(($half_bin_number+$j)) && z=$(($k))
            cat $output_dir/${sample}.$type.${list}.bin_$x.intersect.bed $output_dir/${sample}.$type.${list}.bin_$y.intersect.bed > \
            $output_dir/merge/${sample}.$type.${list}.bin_$z.intersect.bed &
            cat $output_dir/${sample}.$type.bin_$x.txt $output_dir/${sample}.$type.bin_$y.txt > $output_dir/merge/${sample}.$type.bin_$z.txt &
            i=$(($i+1)) && j=$(($j+1)) && k=$(($k-1))
        done
    done
done
wait
    
    merge_output_dir=$output_dir/merge
    [ -f $merge_output_dir/binned.level.txt ] && rm -rf $merge_output_dir/binned.level.txt
    for sample in ${SAMPLE};do
        for mA_list in ${MA_LIST};do
            for j in nucleosome linker;do
                for i in $(seq 1 $half_bin_number);do
                    BASENAME=$(basename "$mA_list")
                    list_tmp=${BASENAME%%.*} && list=${list_tmp##*_}
                    mA=$(less "$merge_output_dir"/"${sample}".$j.${list}.bin_$i.intersect.bed | wc -l)
                    A=$( sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' )
                    C=$( sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' )
                    G=$( sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' )
                    T=$( sed '1d' "$merge_output_dir"/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' )
                    base_number=$( sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $13} END {print sum}' )
                    AT=$( sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' )
                    GC=$( sed '1d' $merge_output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' )
                    FreqSum=$(less $merge_output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $8} END {print sum}')
                    echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${C}\t${G}\t${T}\t${base_number}\t${AT}\t${GC}\t${FreqSum}" >> $merge_output_dir/binned.level.txt
                done
            done 
        done
    done
### EO　merge different samples 

echo "done"