################################################
#File Name: run.optimize.method_2.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 24 Feb 2022 10:14:05 PM CST
################################################

#!/bin/sh 


set -e


[ -d /home/user/data/lit/project/6mA/feature/iNPS/method_2 ] || mkdir -p /home/user/data/lit/project/6mA/feature/iNPS/method_2
[ -d log/iNPS/method_2 ] || mkdir -p log/iNPS/method_2
output_dir=/home/user/data/lit/project/6mA/feature/iNPS/method_2
log_dir=log/iNPS/method_2
ce11_dir=/home/user/data/lit/database/public/genome/ce11/
ce11_chrSize=$ce11_dir/ce11.fasta.size


<<'!'
### nucleosome bed ###

for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    sort -k1,1 -k2,2n /home/user/data/lit/project/6mA/feature/iNPS/method_1/${sample}.nucleosome.bed | bedtools merge -d 50 -c 1 -o count -i - | \
    awk -v OFS='\t' '$4==1{print $1,$2,$3}' > $output_dir/${sample}.nucleosome.bed
done

### nucleosome bed ###


### linker bed ###
for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    sort -k1,1 -k2,2n /home/user/data/lit/project/6mA/feature/iNPS/method_1/${sample}.nucleosome.bed | bedtools merge -d 50 -c 1 -o count -i - | \
    awk -v OFS='\t' '$4==1{print $1,$2-5*5,$2}' > $output_dir/${sample}.linker.left.bed

    sort -k1,1 -k2,2n /home/user/data/lit/project/6mA/feature/iNPS/method_1/${sample}.nucleosome.bed | bedtools merge -d 50 -c 1 -o count -i - | \
    awk -v OFS='\t' '$4==1{print $1,$3,$3+5*5}' > $output_dir/${sample}.linker.right.bed
done
### linker bed ###

### bin ###

for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    echo "bin ${sample}"
    bedtools makewindows -b $output_dir/${sample}.nucleosome.bed  -n 10 -i winnum > \
    $output_dir/${sample}.nucleosome.binned.bed  2> $log_dir/${sample}.nucleosome.bin.log &

    bedtools makewindows -b $output_dir/${sample}.linker.left.bed  -n 5 -i winnum > \
    $output_dir/${sample}.linker.binned.left.bed 2> $log_dir/${sample}.linker.bin.left.log &

    bedtools makewindows -b $output_dir/${sample}.linker.right.bed  -n 5 -i winnum > \
    $output_dir/${sample}.linker.binned.right.bed 2> $log_dir/${sample}.linker.bin.right.log &
done

wait


### bin ###

for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    less $output_dir/${sample}.linker.binned.right.bed | awk -v OFS='\t' '{print $1,$2,$3,$4+5}' > $output_dir/${sample}.linker.binned.right.add5.bed
    cat $output_dir/${sample}.linker.binned.left.bed $output_dir/${sample}.linker.binned.right.add5.bed > $output_dir/${sample}.linker.bed
done

### 6mA level ###

    [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        for j in nucleosome linker;do
            for i in `seq 1 10`;do
                echo "*** process ${sample}.$j.bin_$i"
                less $output_dir/${sample}.$j.binned.bed | \
                        grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
                            $output_dir/${sample}.$j.bin_$i.bed

                bedtools intersect -wb -a $output_dir/${sample}.$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/${sample}.$j.bin_$i.intersect.bed &

                bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/${sample}.$j.bin_$i.bed > \
                $output_dir/${sample}.$j.bin_$i.fasta &
            done
            wait
        done 
    done ) 


    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        for j in nucleosome linker;do
            for i in `seq 1 10`;do
                mA=`less $output_dir/${sample}.$j.bin_$i.intersect.bed | wc -l`
                A=`grep -i -o a $output_dir/${sample}.$j.bin_$i.fasta | wc -l`
                FreqSum=`less $output_dir/${sample}.$j.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`
                echo -e "${sample}.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> $output_dir/binned.level.txt
            done
        done
    done

### 6mA level ###

# 2022-02-24
# nohup bash run.optimize.method_2.sh >>log/iNPS/method_2/run.log 2>&1 &
!

################################

# 2022-02-25 check

# /home/user/data/lit/project/6mA/feature/iNPS/method_2/PC10_LY_MCC-M1.linker.binned.bed: No such file or directory
# line 62 wrong 


### bin ###

for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    cat $output_dir/${sample}.linker.binned.left.bed $output_dir/${sample}.linker.binned.right.add5.bed > $output_dir/${sample}.linker.binned.bed
done

### 6mA level ###


    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        for j in linker;do
            for i in `seq 1 10`;do
                echo "*** process ${sample}.$j.bin_$i"
                less $output_dir/${sample}.$j.binned.bed | \
                        grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
                            $output_dir/${sample}.$j.bin_$i.bed

                bedtools intersect -wb -a $output_dir/${sample}.$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/${sample}.$j.bin_$i.intersect.bed &

                bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/${sample}.$j.bin_$i.bed > \
                $output_dir/${sample}.$j.bin_$i.fasta &
            done
            wait
        done 
    done ) 


    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        for j in linker;do
            for i in `seq 1 10`;do
                mA=`less $output_dir/${sample}.$j.bin_$i.intersect.bed | wc -l`
                A=`grep -i -o a $output_dir/${sample}.$j.bin_$i.fasta | wc -l`
                FreqSum=`less $output_dir/${sample}.$j.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`
                echo -e "${sample}.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> $output_dir/binned.level.txt
            done
        done
    done

### 6mA level ###


# 2022-02-24
# nohup bash run.optimize.method_2.sh >>log/iNPS/method_2/run.log 2>&1 &