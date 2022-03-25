################################################
#File Name: bin/run.C.reinhardti.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 16 Mar 2022 08:07:52 PM CST
################################################

#!/bin/sh 


set -e

###  feature analysis ###

    [ -d /home/user/data2/lit/project/6mA/reproduce/output/feature/iNPS/C.reinhardti ] || mkdir -p /home/user/data2/lit/project/6mA/reproduce/output/feature/iNPS/C.reinhardti
    [ -d log/feature/iNPS/C.reinhardti ] || mkdir -p log/feature/iNPS/C.reinhardti
    output_dir=/home/user/data2/lit/project/6mA/reproduce/output/feature/iNPS/C.reinhardti
    log_dir=log/feature/iNPS/C.reinhardti
    genome_dir=/home/user/data/lit/database/public/genome/C.reinhardtii
    genome_fasta=${genome_dir}/Creinhardtii_281_v5.0.fa
    genome_chrSize=$genome_dir/Creinhardtii_281_v5.0.fa.size
    sort -k1,1 -k2,2n ${genome_chrSize} > $genome_dir/Creinhardtii_281_v5.0.fa.sorted.size
    genome_sorted_chrSize=$genome_dir/Creinhardtii_281_v5.0.fa.sorted.size

    ### nucleosome bed ###
    time (
    for sample in C.reinhardtii;do
        echo "*** generating ${sample}.nucleosome.bed"
        less /home/user/data/lit/project/6mA/reproduce/iNPS_output/${sample}_Gathering.like_bed | sed '/chromosome=/d' | sed '1d' | \
        # nucleosomes are all kept
        # sorted for bedtools complement 
        awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $output_dir/${sample}.nucleosome.bed

        nucleosome_number=`less $output_dir/${sample}.nucleosome.bed | wc -l`
        echo "${sample} nucleosome_number: ${nucleosome_number}"
    done )
    ### nucleosome bed ###
    
    ### linker bed ###
    time (
    for sample in C.reinhardtii;do
        echo "*** generating ${sample}.linker.bed"
        # do not remove the first linker and the last linker
        bedtools complement  -i $output_dir/${sample}.nucleosome.bed -g $genome_sorted_chrSize > \
        $output_dir/${sample}.linker.bed

        linker_number=`less $output_dir/${sample}.linker.bed | wc -l`
        echo "${sample} linker_number: ${linker_number}"
    done )

    ### linker bed ###


    ### 2. binned linker or nucleosome enrichment ###

    ### 2.1 bin ###
    time (
        for sample in C.reinhardtii;do
        echo "bin ${sample}"
            bedtools makewindows -b $output_dir/${sample}.nucleosome.bed  -n 10 -i winnum > \
            $output_dir/${sample}.nucleosome.binned.bed  2> $log_dir/${sample}.nucleosome.bin.log &
            bedtools makewindows -b $output_dir/${sample}.linker.bed  -n 10 -i winnum > \
            $output_dir/${sample}.linker.binned.bed 2> $log_dir/${sample}.linker.bin.log &
        done 
        wait )

    ### 2.1 bin ###

    ### 2.2 6mA level ###
        [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt
        time (
        for sample in C.reinhardtii;do
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


        for sample in C.reinhardtii;do
                for j in nucleosome linker;do
                    for i in `seq 1 10`;do
                        A=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' `
                        C=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' `
                        G=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' `
                        T=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' `
                        AT=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' `
                        GC=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' `

                        echo -e "${sample}\t$j\tbin_$i\t${A}\t${T}\t${G}\t${C}\t${AT}\t${GC}" >> $output_dir/binned.level.txt
                    done
                done
            done
    ### 2.2 6mA level ###

    ### 2. binned linker or nucleosome enrichment ###

### feature analysis ###

# 2022-03-16
# [ -d log/feature/iNPS/C.reinhardti ] || mkdir -p log/feature/iNPS/C.reinhardti
# nohup bash bin/run.C.reinhardti.sh > log/feature/iNPS/C.reinhardti/run.log 2>&1 &