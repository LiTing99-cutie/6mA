################################################
#File Name: run.list_1.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 22 Feb 2022 09:50:45 PM CST
################################################

#!/bin/sh 

set -e

output_dir=/home/user/data/lit/project/6mA/iNPS_output
ce11_dir=/home/user/data/lit/database/public/genome/ce11/


### 2.2 6mA level ###
[ -f binned.level.list_1.txt ] && rm -rf binned.level.list_1.txt
time (
# for i in 1;do
for i in 1 2 3 4 5;do
    for j in nucleosome linker;do
    # for j in nucleosome;do
        for id in 1 2;do
        # for id in 1;do
            echo "*** process PC10_LY_MCC-M${id}.$j.bin_$i"
            # less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.binned.reorg.bed | \
            #         grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
            #             $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed

            bedtools intersect -wb -a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt2.bed | \
            awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.list_1.bed 

            # already get
            # bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed > \
            # $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.fasta

            mA=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.list_1.bed | wc -l`
            A=`grep -i -o a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.fasta | wc -l`
            FreqSum=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.list_1.bed | awk '{sum += $4} END {print sum}'`
            echo -e "PC10_LY_MCC-M${id}.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> binned.level.list_1.txt
        done
    done
done ) 

[ -f binned.level.list_2.txt ] && rm -rf binned.level.list_2.txt
time (
# for i in 1;do
for i in 1 2 3 4 5;do
    for j in nucleosome linker;do
    # for j in nucleosome;do
        for id in 1 2;do
        # for id in 1;do
            echo "*** process PC10_LY_MCC-M${id}.$j.bin_$i"
            # less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.binned.reorg.bed | \
            #         grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
            #             $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed

            bedtools intersect -wb -a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed | \
            awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.list_2.bed 

            # already get
            # bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed > \
            # $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.fasta

            mA=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.list_2.bed | wc -l`
            A=`grep -i -o a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.fasta | wc -l`
            FreqSum=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.list_2.bed | awk '{sum += $4} END {print sum}'`
            echo -e "PC10_LY_MCC-M${id}.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> binned.level.list_2.txt
        done
    done
done ) 


# nohup bash run.list_1.sh >log/feature.list_1_2.log 2>&1 &