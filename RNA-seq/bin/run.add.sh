################################################
#File Name: bin/run.add.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 03 Mar 2022 04:18:47 PM CST
################################################

#!/bin/sh 

for i in `seq 11 12`;do
[ -f data/sample${i}_1.fq.gz ] || ln -s /home/user/data/lit/project/6mA/data/RNA-seq/Clean/RNA_${i}/*_1.fq.gz data/sample${i}_1.fq.gz
[ -f data/sample${i}_2.fq.gz ] || ln -s /home/user/data/lit/project/6mA/data/RNA-seq/Clean/RNA_${i}/*_2.fq.gz data/sample${i}_2.fq.gz
done


### 2. hisat2 mapping ###

    # 2.2 hisat2 mapping 
    hisat2Index=/home/user/data/lit/database/public/genome/ce11/hisat2/ce11/ce11
    fastq_dir=data
    output_dir=/home/user/data/lit/project/6mA/RNA-seq/output/mapping

    for rep in `seq 11 12` ;do
        echo "*** Processing sample${rep}.fastq.gz"
        time (hisat2 -x $hisat2Index -q  \
        -1 $fastq_dir/sample${rep}_1.fq.gz \
        -2 $fastq_dir/sample${rep}_2.fq.gz \
        -S ${output_dir}/sample${rep}.sam \
        --threads 100  1>log/hisat2.sample${rep}.mapping.log 2>&1 
        samtools view -buSh -@ 30 ${output_dir}/sample${rep}.sam > ${output_dir}/sample${rep}.bam)
    done 

### EO 2. hisat2 mapping ###


### 3. sort ###


script_sortBam=/home/user/data2/lit/software/scripts/sortBam.sh 
bam_dir=/home/user/data/lit/project/6mA/RNA-seq/output/mapping

for i in sample;do for j in `seq 11 12`;do echo ${i}${j} ;done ;done |
# head -n 1 | \
while read sample;do
    echo "***`date "+%Y-%m-%d %H:%M:%S"` SAMTOOLS SORT. *** Processing $sample.bam"
    time (bash $script_sortBam \
    -b $bam_dir/$sample.bam \
    -o $bam_dir/$sample.sorted.bam \
    -t 80 \
    -n $sample \
    -ld log)
done 


### 4. mapping qc ###

# 4.1 bam_stat before qc
. /home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh
conda activate RNA-seq
[ -d output/mapping_qc ] || mkdir -p output/mapping_qc
output_dir=output/mapping_qc

for i in sample;do for j in `seq 11 12`;do echo ${i}${j} ;done ;done |
# head -n 1 | \
while read sample;do
    #sorted bam has smaller size
    echo "***`date "+%Y-%m-%d %H:%M:%S"` RSeQC $sample.sorted.bam"
    time (bam_stat.py -i ${bam_dir}/$sample.bam > ${output_dir}/bam_stat_$sample.txt)
done


# 4.2 mark duplicates
script_remDup=/home/user/data2/lit/software/scripts/sambamba_remDup.sh

for i in sample;do for j in `seq 11 12`;do echo ${i}${j} ;done ;done |
# head -n 1 | \
while read sample;do
    echo -e "***`date "+%Y-%m-%d %H:%M:%S"` Mark Duplicates for $sample.sorted.bam"
    time (sambamba markdup -t 10 ${bam_dir}/${sample}.sorted.bam ${bam_dir}/${sample}.markDup.bam
    samtools flagstat -@ 80 ${bam_dir}/${sample}.markDup.bam > ${output_dir}/samtoolsFlagstat.${sample}.markDup.txt)
done
!

### 4.3 remove duplicates and other quality control
. /home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh
conda activate RNA-seq
bam_dir=/home/user/data/lit/project/6mA/RNA-seq/output/mapping
output_dir=output/mapping_qc


for i in sample;do for j in `seq 11 12`;do echo ${i}${j} ;done ;done |
# head -n 1 | \
while read sample;do
    echo "***`date "+%Y-%m-%d %H:%M:%S"` Remove Duplicates and other qc for $sample.markDup.bam"
    time (samtools view -h -f 34 -F 1564 -@ 80  ${bam_dir}/${sample}.markDup.bam > ${bam_dir}/${sample}.clean.1.bam
    samtools view -h -f 18 -F 1580 -@ 80  ${bam_dir}/${sample}.markDup.bam > ${bam_dir}/${sample}.clean.2.bam
    samtools merge -@ 80 -f ${bam_dir}/${sample}.clean.bam ${bam_dir}/${sample}.clean.1.bam ${bam_dir}/${sample}.clean.2.bam 
    samtools flagstat -@ 80 ${bam_dir}/${sample}.clean.bam > ${output_dir}/samtoolsFlagstat.${sample}.clean.txt)
done 


# 4.4 bam_stat after qc

for i in sample;do for j in `seq 11 12`;do echo ${i}${j} ;done ;done |
# head -n 1 | \
while read sample;do
    #sorted bam has smaller size
    echo "***`date "+%Y-%m-%d %H:%M:%S"` RSeQC $sample.clean.bam"
    time (bam_stat.py -i ${bam_dir}/$sample.clean.bam > ${output_dir}/bam_stat_$sample.clean.txt)
done 

### 4. mapping qc ###


### 5. featureCounts ###

[ -d output/featureCounts/add ] || mkdir -p output/featureCounts/add
output_dir=output/featureCounts/add
script_featureCounts=/home/user/data2/lit/software/scripts/featureCounts.sh

echo "*** featureCounts"

# seq 1 for test
SAMPLE=$( for i in sample;do for j in `seq 11 12`;do bam=/home/user/data/lit/project/6mA/RNA-seq/output/mapping/${i}${j}.clean.bam; 
printf "$bam|" ; done ; done )

# Gzipped file is also accepted.
gtf=/home/user/data/lit/database/public/annotation/gtf/ensembl/ensembl.105/Caenorhabditis_elegans.WBcel235.105.gtf.gz

time (bash $script_featureCounts \
-p \
-c \
-g  $gtf  \
-b  $SAMPLE \
-o  $output_dir/featureCounts.txt \
-t  50 \
-n  sample \
-ld log)

### 5. feature counts ###

# nohup bash bin/run.add.sh >> log/run.add.log 2>&1 &



