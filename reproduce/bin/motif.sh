################################################
#File Name: /home/user/data2/lit/project/6mA/reproduce/bin/motif.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 08 Apr 2022 07:54:09 PM CST
################################################

#!/bin/sh 


set -eou pipefail

output_dir=$1
log_dir=$2
center_number=$3
centerORnot=$4
genome_fasta=$5
# need to be sorted first
genome_sorted_chrSize=$6
iNPS_output_dir=$7
SAMPLE=$8

[ -d $output_dir ] || mkdir -p $output_dir
[ -d $log_dir ] || mkdir -p $log_dir

if [ !centerORnot ] ;then
    # nucleosome bed 
    time (
    for sample in $SAMPLE;do
        echo "*** generating ${sample}.nucleosome.bed"
        # clean file no header
        less ${iNPS_output_dir}/${sample}_Gathering.like_bed | \
        awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $output_dir/${sample}.nucleosome.bed
        nucleosome_number=`less $output_dir/${sample}.nucleosome.bed | wc -l`
        echo "${sample} nucleosome_number: ${nucleosome_number}"
    done )

    # linker bed 
    time (
    for sample in $SAMPLE;do
        echo "*** generating ${sample}.linker.bed"
        bedtools complement  -i $output_dir/${sample}.nucleosome.bed -g ${genome_sorted_chrSize} > \
        $output_dir/${sample}.linker.bed

        linker_number=`less $output_dir/${sample}.linker.bed | wc -l`
        echo "${sample} linker_number: ${linker_number}"
    done )

    # binned linker or nucleosome enrichment 

    ## bin
    time (
        for sample in $SAMPLE;do
        echo "bin ${sample}"
            bedtools makewindows -b $output_dir/${sample}.nucleosome.bed  -n 10 -i winnum > \
            $output_dir/${sample}.nucleosome.binned.bed  2> $log_dir/${sample}.nucleosome.bin.log &
            bedtools makewindows -b $output_dir/${sample}.linker.bed  -n 10 -i winnum > \
            $output_dir/${sample}.linker.binned.bed 2> $log_dir/${sample}.linker.bin.log &
        done 
        wait )

    ## 6mA level 
    [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt
    time (
    for sample in $SAMPLE;do
            for j in nucleosome linker;do
                for i in `seq 1 10`;do
                    echo "*** process ${sample}.$j.bin_$i"
                    less $output_dir/${sample}.$j.binned.bed | \
                            grep $i$ > \
                                $output_dir/${sample}.$j.bin_$i.bed

                    bedtools nuc -fi $genome_fasta -bed $output_dir/${sample}.$j.bin_$i.bed -seq > \
                    $output_dir/${sample}.$j.bin_$i.txt &

                done
                wait
            done 
    done ) 


    for sample in $SAMPLE;do
        for j in nucleosome linker;do
            for i in `seq 1 10`;do
                A=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' `
                C=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' `
                G=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' `
                T=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' `
                base_number=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $13} END {print sum}' `
                echo -e "${sample}\t$j\tbin_$i\t${A}\t${C}\t${G}\t${T}\t${base_number}" >> $output_dir/binned.level.txt
            done
        done
    done
fi

if [ centerORnot ] ;then
    echo "happy"
fi
