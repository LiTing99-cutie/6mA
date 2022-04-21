################################################
#File Name: bin/run.motif.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sun 13 Mar 2022 08:54:16 PM CST
################################################

#!/bin/sh 

<<'!'
# 1. get nucleosome or linker fasta
[ -d /home/user/data/lit/project/6mA/feature/iNPS/motif ] || mkdir -p /home/user/data/lit/project/6mA/feature/iNPS/motif
[ -d log/iNPS/motif ] || mkdir -p log/iNPS/motif
input_dir=/home/user/data/lit/project/6mA/feature/iNPS/method_2
output_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif
ce11_dir=/home/user/data/lit/database/public/genome/ce11/

for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    echo "get nucleosome fasta"
    bedtools getfasta -fi $ce11_dir/ce11.fa -bed $input_dir/${sample}.nucleosome.bed > \
                $input_dir/${sample}.nucleosome.fasta &
    echo "get linker fasta"
    cat $input_dir/${sample}.linker.left.bed $input_dir/${sample}.linker.right.bed > $input_dir/${sample}.linker.bed
    bedtools getfasta -fi $ce11_dir/ce11.fa -bed $input_dir/${sample}.linker.bed > \
                $input_dir/${sample}.linker.fasta &
done 

# nohup bash bin/run.motif.sh > log/iNPS/motif/getfasta.log 2>&1 &
!

# install seqkit ( locate function )

# conda install -c bioconda csvtk
# conda install -c bioconda seqkit

# -i -ignore-case
# -p pattern
# -j threads

<<'!'
# 2. extract motif
# 3-base motif
input_dir=/home/user/data/lit/project/6mA/feature/iNPS/method_2
output_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif

for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
for j in nucleosome linker;do
    echo "extract motif for ${sample}.${j}"
    time ( cat $input_dir/${sample}.${j}.fasta | \
    # modify -p here
    seqkit locate -i -d -p NAN -j 20 | csvtk pretty -t | sed '1,2d' | awk -F '[\t:]' -v OFS='\t' '{print $1,$2}' | \
    awk -v OFS='\t' '{print $1,($6+$7)/2-1,($6+$7)/2,$8}' > ${output_dir}/${sample}.${j}.motif.bed )
done
done


# 3. stat motif feature

[ -f $output_dir/level.txt ] && rm -rf $output_dir/level.txt
time (
for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
for j in nucleosome linker; do

        list=${mA_list%.*} && list=${mA_list##*_}   
        echo "stat for ${sample}.${j}.${list} "

        bedtools intersect -wa -wb -a ${output_dir}/${sample}.${j}.motif.bed -b ${mA_list} | \
        awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' > $output_dir/${sample}.${j}.${list}.intersect.bed 

        awk '{print $4}' $output_dir/${sample}.${j}.bed | sort | uniq > $output_dir/${sample}.${j}.motif.txt

        cat $output_dir/${sample}.${j}.motif.txt | \
        while read motif;do | \
            awk '{if($4=="'$motif'") print$0}' $output_dir/${sample}.${j}.bed > $output_dir/${sample}.${j}.${motif}.bed
            awk '{if($4=="'$motif'") print$0}' $output_dir/${sample}.${j}.${list}.intersect.bed > $output_dir/${sample}.${j}.${list}.${motif}.intersect.bed
            FreqSum=`less $output_dir/${sample}.${j}.${list}.${motif}.intersect.bed | awk '{sum += $5} END {print sum}'`
            mA=`less $output_dir/${sample}.${j}.${list}.${motif}.intersect.bed | wc -l`
            A=`less ${output_dir}/${sample}.${j}.${motif}.bed | wc -l`
            
            echo -e "${sample}\t${j}\t${motif}\t${mA}\t${A}\t${FreqSum}\t${list}" >> $output_dir/level.txt
        done
done
done
done )

# nohup bash bin/run.motif.sh > log/iNPS/motif/part.2-3.log 2>&1 &
!

# line 72 wrong

# 3. stat motif feature

input_dir=/home/user/data/lit/project/6mA/feature/iNPS/method_2
output_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif

[ -f $output_dir/level.txt ] && rm -rf $output_dir/level.txt
time (
for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
for j in nucleosome linker; do

        list=${mA_list%.*} && list=${mA_list##*_}   
        echo "stat for ${sample}.${j}.${list} "

        bedtools intersect -wa -wb -a ${output_dir}/${sample}.${j}.motif.bed -b ${mA_list} | \
        awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' > $output_dir/${sample}.${j}.${list}.motif.intersect.bed 

        awk '{print $4}' $output_dir/${sample}.${j}.motif.bed | sort | uniq > $output_dir/${sample}.${j}.motif.txt

        cat $output_dir/${sample}.${j}.motif.txt | \
        while read motif;do 
            awk '{if($4=="'$motif'") print$0}' $output_dir/${sample}.${j}.motif.bed > $output_dir/${sample}.${j}.${motif}.bed
            awk '{if($4=="'$motif'") print$0}' $output_dir/${sample}.${j}.${list}.motif.intersect.bed > $output_dir/${sample}.${j}.${list}.${motif}.intersect.bed
            FreqSum=`less $output_dir/${sample}.${j}.${list}.${motif}.intersect.bed | awk '{sum += $5} END {print sum}'`
            mA=`less $output_dir/${sample}.${j}.${list}.${motif}.intersect.bed | wc -l`
            A=`less ${output_dir}/${sample}.${j}.${motif}.bed | wc -l`
            
            echo -e "${sample}\t${j}\t${motif}\t${mA}\t${A}\t${FreqSum}\t${list}" >> $output_dir/level.txt
        done
done
done
done )

# nohup bash bin/run.motif.sh > log/iNPS/motif/part.3.log 2>&1 &