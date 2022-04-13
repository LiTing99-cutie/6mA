################################################
#File Name: bin/mnase.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 08 Mar 2022 06:46:49 PM CST
################################################

#!/bin/sh 



# set -e

<<'!'
# 1. merge data and create soft links
Mnase_data_dir=/home/user/data2/lit/DATA/project/6mA/reproduce/rawdata/Mnase
input_dir=input/rawdata/Mnase

for i in 1 2;do
    zcat $Mnase_data_dir/SRR1994849_GSM1666100_MNase-seq_Chlamydomonas_reinhardtii_MNase-Seq_$i.fastq.gz \
    $Mnase_data_dir/SRR1994850_GSM1666100_MNase-seq_Chlamydomonas_reinhardtii_MNase-Seq_$i.fastq.gz | \
    gzip > $Mnase_data_dir/C.reinhardtii_$i.fastq.gz
    [ -f $input_dir/T.thermophila_$i.fastq.gz ] || ln -s $Mnase_data_dir/SRR5337752_GSM2534785_MNase-seq_Tetrahymena_thermophila_MNase-Seq_$i.fastq.gz \
    $input_dir/T.thermophila_$i.fastq.gz
    [ -f $input_dir/C.reinhardtii_$i.fastq.gz ] || ln -s $Mnase_data_dir/C.reinhardtii_$i.fastq.gz \
    $input_dir/C.reinhardtii_$i.fastq.gz
done

# 2. rawdata quality control (fastqc) 
    . /home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh
    conda activate RNA-seq

    [ -d output/fastqc_output/ ] || mkdir -p output/fastqc_output/
    [ -d output/multiqc_output/ ] || mkdir -p output/multiqc_output/

    fastqc_output_dir=`pwd`/output/fastqc_output/
    multiqc_output_dir=`pwd`/output/multiqc_output/

    pushd $input_dir
    ls | \
    # head -n1 | \
    while read fastq;do
    fastqc -t 50 -o $fastqc_output_dir $fastq
    done
    popd

    multiqc -o $multiqc_output_dir $fastqc_output_dir

    # nohup bash bin/mnase.sh > log/mnase.fastqc.log 2>&1 &
!

<<'!'
# 3. build bowtie indexes (bowtie2)
    [ -d log ] || mkdir log  
    [ -d /home/user/data/lit/project/6mA/reproduce/mapping ] || mkdir /home/user/data/lit/project/6mA/reproduce/mapping
    [ -d output/qc ] || mkdir output/qc

    # bowtie2-build accept gzip files
    for spe in C.reinhardtii T.thermophila;do
        echo "***`date "+%Y-%m-%d %H:%M:%S"` bowtie2-build $spe"
        [ -d /home/user/data2/lit/DATA/database/public/genome/$spe/index/bowtie2/ ] || \
        mkdir -p /home/user/data2/lit/DATA/database/public/genome/$spe/index/bowtie2/
        genome_fasta=/home/user/data/lit/database/public/genome/$spe/*fa*
        time (bowtie2-build --threads 20 $genome_fasta \
        /home/user/data2/lit/DATA/database/public/genome/$spe/index/bowtie2/$spe) 
    done

# 4. mapping (bowtie2)
    input_dir=input/rawdata/Mnase
    [ -d /home/user/data2/lit/DATA/project/6mA/reproduce/mapping ] || mkdir /home/user/data2/lit/DATA/project/6mA/reproduce/mapping
    mapping_dir=/home/user/data2/lit/DATA/project/6mA/reproduce/mapping
    for sample in C.reinhardtii T.thermophila;do
        echo "***`date "+%Y-%m-%d %H:%M:%S"` mapping $sample"
        fq_1=$input_dir/${sample}_1.fastq.gz
        fq_2=$input_dir/${sample}_2.fastq.gz
        spe=$sample
        time (bowtie2 -p 20 -x /home/user/data2/lit/DATA/database/public/genome/$spe/index/bowtie2/$spe \
        -1 $fq_1 -2 $fq_2 -S $mapping_dir/${sample}.sam) 
    done
!

<<'!'
### 5. mapping qc ###

# 5.0 sam to bam and sort
script_sortBam=/home/user/data2/lit/software/scripts/sortBam.sh 
bam_dir=/home/user/data2/lit/DATA/project/6mA/reproduce/mapping

for i in C.reinhardtii T.thermophila;do echo ${i};done |
# head -n 1 | \
while read sample;do
    echo "***`date "+%Y-%m-%d %H:%M:%S"` SAMTOOLS SORT. *** Processing $sample.bam"
    samtools view -buSh ${bam_dir}/$sample.sam > ${bam_dir}/$sample.bam
    time (bash $script_sortBam \
    -b $bam_dir/$sample.bam \
    -o $bam_dir/$sample.sorted.bam \
    -t 80 \
    -n $sample \
    -ld log)
done 

# 5.1 bam_stat before qc
. /home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh
conda activate RNA-seq
[ -d output/mapping_qc ] || mkdir -p output/mapping_qc
output_dir=output/mapping_qc

for i in C.reinhardtii T.thermophila;do echo ${i};done |
# head -n 1 | \
while read sample;do
    #sorted bam has smaller size
    echo "***`date "+%Y-%m-%d %H:%M:%S"` RSeQC $sample.sorted.bam"
    time (bam_stat.py -i ${bam_dir}/$sample.sorted.bam > ${output_dir}/bam_stat_$sample.txt) 
done

# 5.2 mark duplicates
script_remDup=/home/user/data2/lit/software/scripts/sambamba_remDup.sh

for i in C.reinhardtii T.thermophila;do echo ${i};done |
# head -n 1 | \
while read sample;do
    echo -e "***`date "+%Y-%m-%d %H:%M:%S"` Mark Duplicates for $sample.sorted.bam"
    time (sambamba markdup -t 10 ${bam_dir}/${sample}.sorted.bam ${bam_dir}/${sample}.markDup.bam
    samtools flagstat -@ 80 ${bam_dir}/${sample}.markDup.bam > ${output_dir}/samtoolsFlagstat.${sample}.markDup.txt) 
done


### 5.3 remove duplicates and other quality control

for i in C.reinhardtii T.thermophila;do echo ${i};done |
# head -n 1 | \
while read sample;do
    echo "***`date "+%Y-%m-%d %H:%M:%S"` Remove Duplicates and other qc for $sample.markDup.bam"
    time (samtools view -h -f 34 -F 1564 -@ 80  ${bam_dir}/${sample}.markDup.bam > ${bam_dir}/${sample}.clean.1.bam
    samtools view -h -f 18 -F 1580 -@ 80  ${bam_dir}/${sample}.markDup.bam > ${bam_dir}/${sample}.clean.2.bam
    samtools merge -@ 80 -f ${bam_dir}/${sample}.clean.bam ${bam_dir}/${sample}.clean.1.bam ${bam_dir}/${sample}.clean.2.bam 
    samtools flagstat -@ 80 ${bam_dir}/${sample}.clean.bam > ${output_dir}/samtoolsFlagstat.${sample}.clean.txt)
done 


# 5.4 bam_stat after qc

for i in C.reinhardtii T.thermophila;do echo ${i};done |
# head -n 1 | \
while read sample;do
    #sorted bam has smaller size
    echo "***`date "+%Y-%m-%d %H:%M:%S"` RSeQC $sample.clean.bam"
    time (bam_stat.py -i ${bam_dir}/$sample.clean.bam > ${output_dir}/bam_stat_$sample.clean.txt)
done


# 6. iNPS 
    [ -d /home/user/data/lit/project/6mA/reproduce/iNPS_output ] || mkdir -p /home/user/data/lit/project/6mA/reproduce/iNPS_output
    output_dir=/home/user/data/lit/project/6mA/reproduce/iNPS_output
    mapping_dir=/home/user/data/lit/project/6mA/reproduce/mapping
     
    for i in C.reinhardtii T.thermophila;do echo ${i};done |
    # head -n 1 | \
    while read sample;do
        echo "***`date "+%Y-%m-%d %H:%M:%S"` iNPS $sample"
        samtools sort -n -O BAM -@ 80 $mapping_dir/${sample}.clean.bam > $mapping_dir/${sample}.clean.sorted.bam
        bedtools bamtobed -bedpe -i $mapping_dir/${sample}.clean.sorted.bam > $mapping_dir/${sample}.bedpe
        awk  -v OFS='\t' '{print $1,$2,$6}'  $mapping_dir/${sample}.bedpe > $mapping_dir/${sample}.bedpe.bed
        iNPS=/home/user/data2/lit/software/bin/iNPS_V1.2.2.py
        time (python3 $iNPS --s_p=p -i $mapping_dir/${sample}.bedpe.bed -o $output_dir/${sample})
    done

    # nohup bash bin/mnase.sh > log/run.mnase.3-6.log 2>&1 &
    # nohup bash bin/mnase.sh > log/run.mnase.5-6.log 2>&1 &


multiqc -f  -o output/multiqc_output output/mapping_qc/*clean.txt -n mapping_qc
!

<<'!'
# 6. iNPS 
    [ -d /home/user/data/lit/project/6mA/reproduce/iNPS_output ] || mkdir -p /home/user/data/lit/project/6mA/reproduce/iNPS_output
    output_dir=/home/user/data/lit/project/6mA/reproduce/iNPS_output
    mapping_dir=/home/user/data/lit/project/6mA/reproduce/mapping
     
    for i in C.reinhardtii T.thermophila;do echo ${i};done |
    # head -n 1 | \
    while read sample;do
        samtools sort -n -O BAM -@ 80 $mapping_dir/${sample}.clean.bam > $mapping_dir/${sample}.clean.sorted.bam
        bedtools bamtobed -bedpe -i $mapping_dir/${sample}.clean.sorted.bam > $mapping_dir/${sample}.bedpe
        awk  -v OFS='\t' '{print $1,$2,$6}'  $mapping_dir/${sample}.bedpe > $mapping_dir/${sample}.bedpe.bed
        iNPS=/home/user/data2/lit/software/bin/iNPS_V1.2.2.py
        echo "***`date "+%Y-%m-%d %H:%M:%S"` iNPS $sample"
        # iNPS log file too large
        time ( python3 $iNPS --s_p=p -i $mapping_dir/${sample}.bedpe.bed -o $output_dir/${sample} )
    done

    # nohup bash bin/mnase.sh > log/run.mnase.6.log 2>&1 &
    # nohup bash bin/mnase.sh 2> log/run.mnase.6.error.log 1> /dev/null &
!

<<'!'
# 6. iNPS
    [ -d /home/user/data/lit/project/6mA/reproduce/iNPS_output ] || mkdir -p /home/user/data/lit/project/6mA/reproduce/iNPS_output
    output_dir=/home/user/data/lit/project/6mA/reproduce/iNPS_output
    mapping_dir=/home/user/data/lit/project/6mA/reproduce/mapping

    for i in T.thermophila;do echo ${i};done |
    # head -n 1 | \
    while read sample;do
        # samtools sort -n -O BAM -@ 80 $mapping_dir/${sample}.clean.bam > $mapping_dir/${sample}.clean.sorted.bam
        # bedtools bamtobed -bedpe -i $mapping_dir/${sample}.clean.sorted.bam > $mapping_dir/${sample}.bedpe
        # awk  -v OFS='\t' '{print $1,$2,$6}'  $mapping_dir/${sample}.bedpe > $mapping_dir/${sample}.bedpe.bed
		awk '{print $1}' $mapping_dir/${sample}.bedpe.bed | sort | uniq > $mapping_dir/${sample}.scaffold.txt
		cat $mapping_dir/${sample}.scaffold.txt| \
        # head -n1 | \
		while read scaffold;do
		awk '{if($1=="'${scaffold}'") print $0}' $mapping_dir/${sample}.bedpe.bed > $mapping_dir/${sample}.bedpe.${scaffold}.bed
		iNPS=/home/user/data2/lit/software/bin/iNPS_V1.2.2.py
         
        # iNPS log file too large
        python3 $iNPS --s_p=p -i $mapping_dir/${sample}.bedpe.${scaffold}.bed -o $output_dir/${sample}.${scaffold}
        echo "***`date "+%Y-%m-%d %H:%M:%S"` iNPS ${sample}.${scaffold} $?" >> log/run.mnase.T.log
		done
	done

	# nohup bash bin/mnase.sh 2> log/run.mnase.T.error.log 1> /dev/null &
!

<<'!'
# 6. iNPS 
    # annotate set -e 
    [ -d /home/user/data/lit/project/6mA/reproduce/iNPS_output ] || mkdir -p /home/user/data/lit/project/6mA/reproduce/iNPS_output
    output_dir=/home/user/data/lit/project/6mA/reproduce/iNPS_output
    mapping_dir=/home/user/data/lit/project/6mA/reproduce/mapping
     
    for i in T.thermophila;do echo ${i};done |
    # head -n 1 | \
    while read sample;do
        iNPS=/home/user/data2/lit/software/bin/iNPS_V1.2.2.py
        echo "***`date "+%Y-%m-%d %H:%M:%S"` iNPS $sample"
        # iNPS log file too large
        time ( python3 $iNPS --s_p=p -i $mapping_dir/${sample}.bedpe.bed -o $output_dir/${sample} )
    done

    # nohup bash bin/mnase.sh 2> log/run.mnase.6.error.1.log 1> /dev/null &
!

# 6. iNPS
    # annotate set -e
    [ -d /home/user/data/lit/project/6mA/reproduce/iNPS_output ] || mkdir -p /home/user/data/lit/project/6mA/reproduce/iNPS_output
    output_dir=/home/user/data/lit/project/6mA/reproduce/iNPS_output
    mapping_dir=/home/user/data/lit/project/6mA/reproduce/mapping

    for i in T.thermophila;do echo ${i};done |
    # head -n 1 | \
    while read sample;do
        # samtools sort -n -O BAM -@ 80 $mapping_dir/${sample}.clean.bam > $mapping_dir/${sample}.clean.sorted.bam
        # bedtools bamtobed -bedpe -i $mapping_dir/${sample}.clean.sorted.bam > $mapping_dir/${sample}.bedpe
        # awk  -v OFS='\t' '{print $1,$2,$6}'  $mapping_dir/${sample}.bedpe > $mapping_dir/${sample}.bedpe.bed
		awk '{print $1}' $mapping_dir/${sample}.bedpe.bed | sort | uniq > $mapping_dir/${sample}.scaffold.txt
		cat $mapping_dir/${sample}.scaffold.txt| \
        # head -n1 | \
		while read scaffold;do
		awk '{if($1=="'${scaffold}'") print $0}' $mapping_dir/${sample}.bedpe.bed > $mapping_dir/${sample}.bedpe.${scaffold}.bed
		iNPS=/home/user/data2/lit/software/bin/iNPS_V1.2.2.py
         
        # iNPS log file too large
        python3 $iNPS --s_p=p -i $mapping_dir/${sample}.bedpe.${scaffold}.bed -o $output_dir/${sample}.${scaffold}
        echo "***`date "+%Y-%m-%d %H:%M:%S"` iNPS ${sample}.${scaffold} $?" >> log/run.mnase.T.log
		done
	done

	# nohup bash bin/mnase.sh 2> log/run.mnase.T.error.2.log 1> /dev/null &


# less /home/user/data/lit/project/6mA/reproduce/mapping/T.thermophila.scaffold.txt | wc -l
# 1142

# T.thermophila.scf_8253829_Gathering.like_bed

# ls /home/user/data/lit/project/6mA/reproduce/iNPS_output/T.thermophila*_Gathering.like_bed | wc -l
# 1130

cat /home/user/data/lit/project/6mA/reproduce/iNPS_output/T.thermophila*_Gathering.like_bed | \
sed '/chromosome=/d;/Chromosome/d' > \
/home/user/data/lit/project/6mA/reproduce/iNPS_output/T.thermophila_Gathering.like_bed

# sed '/Chromosome/d'|less