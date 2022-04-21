################################################
#File Name: run.relinker.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 23 Feb 2022 12:08:27 PM CST
################################################

#!/bin/sh 

set -e 

output_dir=/home/user/data/lit/project/6mA/iNPS_output
ce11_dir=/home/user/data/lit/database/public/genome/ce11/
ce11_chrSize=$ce11_dir/ce11.fasta.size
# awk -v OFS="\t" '{print $1,0,$2}' $ce11_chrSize > $ce11_dir/ce11.fasta.size.bed

### nucleosome bed ###

    # already done

### nucleosome bed ###


<<'!'
### linker bed ###
    for i in 1 2;do
        sort -k1,1 -k2,2n PC10_LY_MCC-M${i}_Gathering.nucleosome.bed | bedtools merge -d 50 -c 1 -o count -i - | \
        awk -v OFS='\t' '$4==1{print $1,$2-5*5,$2,"left"}' > PC10_LY_MCC-M${i}_Gathering.linker.left.bed

        sort -k1,1 -k2,2n PC10_LY_MCC-M${i}_Gathering.nucleosome.bed | bedtools merge -d 50 -c 1 -o count -i - | \
        awk -v OFS='\t' '$4==1{print $1,$3,$3+5*5,"right"}' > PC10_LY_MCC-M${i}_Gathering.linker.right.bed
    done
### linker bed ###
!


<<'!'
### bin linker ( nucleosome already bin ) ###
    for i in 1 2;do
        cat PC10_LY_MCC-M${i}_Gathering.linker.left.bed PC10_LY_MCC-M${i}_Gathering.linker.right.bed > PC10_LY_MCC-M${i}_Gathering.linker.n.bed

        bedtools makewindows -b PC10_LY_MCC-M${i}_Gathering.linker.n.bed -n 5 -i srcwinnum > \
        $output_dir/PC10_LY_MCC-M${i}_Gathering.linker.n.binned.bed
    done
### bin linker ( nucleosome already bin ) ###
!

### 6mA level ###

<<'!'
### no multithreading ###
    [ -f binned.level.n.txt ] && rm -rf binned.level.n.txt

    for i in $(for i in `seq 1 5`;do echo left_$i;echo right_$i;done);do
        for j in linker;do
            for id in 1 2;do
                echo "*** process PC10_LY_MCC-M${id}.$j.bin_$i"
                less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.n.binned.bed | \
                        grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
                            $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed
                bedtools intersect -wb -a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.bed 
                bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed > \
                $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.fasta
                mA=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.bed | wc -l`
                A=`grep -i -o a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.fasta | wc -l`
                FreqSum=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`
                echo -e "PC10_LY_MCC-M${id}.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> binned.level.n.txt
            done
        done
    done
### no multithreading ###
!

# 2022-02-24
### Manual multithreading ###
    [ -f binned.level.ten_bins.txt ] && rm -rf binned.level.ten_bins.txt

    time (for id in 1 2;do
        for j in nucleosome;do
            for i in `seq 1 10`;do
                echo "*** process PC10_LY_MCC-M${id}.$j.bin_$i"
                less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.binned.bed | \
                        grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
                            $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.ten_bins.bed

                bedtools intersect -wb -a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.ten_bins.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.ten_bins.bed &

                bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.ten_bins.bed > \
                $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.ten_bins.fasta &
            done
            wait
        done
    done )

    for id in 1 2;do
        for j in nucleosome;do
            for i in `seq 1 10`;do
                mA=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.ten_bins.bed | wc -l`
                A=`grep -i -o a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.ten_bins.fasta | wc -l`
                FreqSum=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.ten_bins.bed | awk '{sum += $4} END {print sum}'`
                echo -e "PC10_LY_MCC-M${id}.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> binned.level.ten_bins.txt
            done
        done
    done

### Manual multithreading ###

### 6mA level ###

# 2022-02-23
# nohup bash run.relinker.sh >log/feature.relinker.log 2>&1 &

# 2022-02-24
# nohup bash run.relinker.sh >>log/feature.relinker.log 2>&1 &
# nohup bash run.relinker.sh >>log/feature.relinker.log 2>&1 &