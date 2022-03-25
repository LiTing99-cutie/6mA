################################################
#File Name: /home/user/data2/lit/project/6mA/run.danpos.n.1_1.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 15 Mar 2022 05:08:30 PM CST
################################################

#!/bin/sh 


set -e

###  feature analysis ###

    [ -d /home/user/data/lit/project/6mA/feature/danpos/normed/method_1_1 ] || mkdir -p /home/user/data/lit/project/6mA/feature/danpos/normed/method_1_1
    [ -d log/danpos/normed/method_1_1 ] || mkdir -p log/danpos/normed/method_1_1
    output_dir=/home/user/data/lit/project/6mA/feature/danpos/normed/method_1_1
    log_dir=log/danpos/normed/method_1_1
    ce11_dir=/home/user/data/lit/database/public/genome/ce11/pengq/
    ce11_chrSize=$ce11_dir/Caenorhabditis_elegans.WBcel235.dna.fa.size

    ### nucleosome bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.nucleosome.bed"
        sed '1d' danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls | \
        # nucleosomes are all kept
        # sorted for bedtools complement 
        awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $output_dir/${sample}.nucleosome.bed

        nucleosome_number=`less $output_dir/${sample}.nucleosome.bed | wc -l`
        echo "${sample} nucleosome_number: ${nucleosome_number}"
    done )
    ### nucleosome bed ###
    
    ### linker bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.linker.bed"
        # do not remove the first linker and the last linker
        bedtools complement  -i $output_dir/${sample}.nucleosome.bed -g $ce11_chrSize > \
        $output_dir/${sample}.linker.bed

        linker_number=`less $output_dir/${sample}.linker.bed | wc -l`
        echo "${sample} linker_number: ${linker_number}"
    done )

    ### linker bed ###


    ### 2. binned linker or nucleosome enrichment ###

    ### 2.1 bin ###
    time (
        for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
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
        for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
            for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
                for j in nucleosome linker;do
                    for i in `seq 1 10`;do
                        list_tmp=${mA_list%.*} && list=${list_tmp##*_} 
                        echo "*** process ${sample}.$j.${list}.bin_$i"
                        less $output_dir/${sample}.$j.binned.bed | \
                                grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
                                    $output_dir/${sample}.$j.bin_$i.bed
                                
                        less $output_dir/${sample}.$j.bin_$i.bed | awk -v OFS='\t' '{if($1=="MtDNA") print "chrM",$2,$3 ; else print "chr"$1,$2,$3}' > \
                                    $output_dir/${sample}.$j.bin_$i.chr.bed

                        bedtools intersect -wb -a $output_dir/${sample}.$j.bin_$i.chr.bed -b ${mA_list} | \
                        awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed &

                        bedtools nuc -fi $ce11_dir/Caenorhabditis_elegans.WBcel235.dna.fa -bed $output_dir/${sample}.$j.bin_$i.bed > \
                        $output_dir/${sample}.$j.bin_$i.txt &

                    done
                    wait
                done 
            done
        done ) 


        for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
            for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
                for j in nucleosome linker;do
                    for i in `seq 1 10`;do
                        list_tmp=${mA_list%.*} && list=${list_tmp##*_}
                        mA=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | wc -l`
                        A=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum}' `
                        T=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' `
                        AT=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $4} END {print sum/NR}' `
                        GC=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' `
                        FreqSum=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`
                        echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${T}\t${FreqSum}\t${AT}\t${GC}" >> $output_dir/binned.level.txt
                    done
                done
            done
        done
    ### 2.2 6mA level ###

    ### 2. binned linker or nucleosome enrichment ###

### feature analysis ###

# 2022-03-16
# [ -d log/danpos/normed/method_1_1 ] || mkdir -p log/danpos/normed/method_1_1
# nohup bash run.danpos.n.1_1.sh > log/danpos/normed/method_1_1/run.log 2>&1 &

<<!
# 2022-03-16 difference between two annotation file
    bedtools getfasta -fi /home/user/data/lit/database/public/genome/ce11/ce11.fa -bed sig_gt_2_fc0_cov10_for_lt_wt1.bed > \
                        wt1.fasta &
    less wt1.fasta | grep -v '^>' | sort | uniq -c  

    #     710 a
    #    6196 A
    #     669 t
    #    5516 T

    less sig_gt_2_fc0_cov10_for_lt_wt1.bed | sed 's/chr//g;s/M/MtDNA/g'| cut -f 1-3 > wt1.chr.rm.bed
    bedtools getfasta -fi /home/user/data/lit/database/public/genome/ce11/pengq/Caenorhabditis_elegans.WBcel235.dna.fa -bed wt1.chr.rm.bed > \
                        wt1.chr.rm.fasta &

    less wt1.chr.rm.fasta | grep -v '^>' | sort | uniq -c 

    #    6906 A
    #    6185 T
!
