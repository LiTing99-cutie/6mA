################################################
#File Name: bin/run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 25 Feb 2022 09:24:31 PM CST
################################################

#!/bin/sh 

<<'!'
### 1. fastqc ###
    for i in `seq 1 10`;do
    [ -f data/sample${i}_1.fq.gz ] || ln -s /home/user/data/lit/project/6mA/data/RNA-seq/Clean/RNA_${i}/*_1.fq.gz data/sample${i}_1.fq.gz
    [ -f data/sample${i}_2.fq.gz ] || ln -s /home/user/data/lit/project/6mA/data/RNA-seq/Clean/RNA_${i}/*_2.fq.gz data/sample${i}_2.fq.gz
    done


    . /home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh
    conda activate RNA-seq

    [ -d output/fastqc_output/ ] || mkdir -p output/fastqc_output/
    [ -d output/multiqc_output/ ] || mkdir -p output/multiqc_output/

    fastqc_output_dir=`pwd`/output/fastqc_output/
    multiqc_output_dir=`pwd`/output/multiqc_output/

    pushd data
    ls | \
    # head -n1 | \
    while read fastq;do
    fastqc -t 50 -o $fastqc_output_dir $fastq
    done
    popd

    multiqc -o $multiqc_output_dir $fastqc_output_dir

    # nohup bash bin/run.sh > log/fastqc.log 2>&1 &
### EO 1. fastqc ###
!

<<'!'
### 2. hisat2 mapping ###

    # 2.1 Build hisat2 index for panTro5 

    # for spe in ce11;do
    #     data_dir=/home/user/data/lit/database/public/genome/${spe}
    #     [ -d $data_dir/hisat2/${spe} ]  ||  mkdir -p $data_dir/hisat2/${spe}
    #     fasta=$data_dir/${spe}.fa
    #     time (hisat2-build $fasta $data_dir/hisat2/${spe}/${spe} 1>log/hisat2-build.${spe}.log 2>&1)
    # done

    # 2.2 hisat2 mapping 
    hisat2Index=/home/user/data/lit/database/public/genome/ce11/hisat2/ce11/ce11
    fastq_dir=data
    [ -d output/mapping/ ] || mkdir -p output/mapping/
    output_dir=output/mapping/

    for rep in `seq 1 10` ;do
        echo "*** Processing sample${rep}.fastq.gz"
        time (hisat2 -x $hisat2Index -q  \
        -1 $fastq_dir/sample${rep}_1.fq.gz \
        -2 $fastq_dir/sample${rep}_2.fq.gz \
        -S ${output_dir}/sample${rep}.sam \
        --threads 100  1>log/hisat2.sample${rep}.mapping.log 2>&1 
        samtools view -buSh -@ 30 ${output_dir}/sample${rep}.sam > ${output_dir}/sample${rep}.bam)
    done 

    # nohup bash bin/run.sh > log/hisat2.log 2>&1 &

    # nohup bash bin/run.sh >> log/hisat2.log 2>&1 &

    # nohup cp -r output/mapping/* /home/user/data/lit/project/6mA/RNA-seq/output/mapping &    

    # rm -rf output/mapping/

### EO 2. hisat2 mapping ###
!

### 3. sort ###


script_sortBam=/home/user/data2/lit/software/scripts/sortBam.sh 
bam_dir=/home/user/data/lit/project/6mA/RNA-seq/output/mapping

for i in sample;do for j in `seq 1 10`;do echo ${i}${j} ;done ;done |
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


<<'!'
### ref ###
    #### 1. Sort bam ; Remove Duplicates ; library size
        script_sortBam=/home/user/data2/lit/software/scripts/sortBam.sh 
        script_remDup=/home/user/data2/lit/software/scripts/sambamba_remDup.sh

        # [ -f log/hrm_developmental/sortBam.picard_remDup.libSize.log ] && rm -rf log/hrm_developmental/sortBam.picard_remDup.libSize.log
        [ -f log/hrm_developmental/sortBam.picard_remDup.libSize.1.log ] && rm -rf log/hrm_developmental/sortBam.picard_remDup.libSize.1.log
        # [ -f $analysis_dir/rpkm/hrm_developmental/mappedReadCounts.txt  ] && rm -f $analysis_dir/rpkm/hrm_developmental/mappedReadCounts.txt
        [ -f $analysis_dir/rpkm/hrm_developmental/mappedReadCounts.1.txt  ] && rm -f $analysis_dir/rpkm/hrm_developmental/mappedReadCounts.1.txt

        # for spe in human;do
        # for spe in human rhesus mouse;do
        for spe in rhesus mouse;do
            bamLst=/home/user/data/uplee/data/dataset/paper/2019-nature-moreira/bam/bam/$spe/*sorted.bam
            ls $bamLst | 
            # head -n1 | 
            while read bam;do
                sample=$( echo $bam | sed -E 's@^.*/@@g;s@.sorted.bam@@g' )
                species=$( echo $bam | awk -F'.' '{print tolower($3)}' | sed 's@macaque@rhesus@g' )
                    
                echo "***`date "+%Y-%m-%d %H:%M:%S"` SAMTOOLS SORT. *** Processing $sample"
                bash $script_sortBam \
                -b $mappingDir/$species/$sample.bam \
                -o $analysis_dir/mapping/hrm_developmental/$species/$sample.lsorted.bam \
                -t 80 \
                -n $sample \
                -ld log
                    
                echo "***`date "+%Y-%m-%d %H:%M:%S"` PICARD REMOVE DUPLICATES. Processing $sample"
                bash $script_remDup \
                -i $analysis_dir/mapping/hrm_developmental/$species/$sample.lsorted.bam  \
                -md $analysis_dir/mapping/hrm_developmental/$species/$sample.markDup.bam  \
                -rd $analysis_dir/mapping/hrm_developmental/$species/$sample.remDup.bam \
                -n $sample \
                -ld log

                echo "***`date "+%Y-%m-%d %H:%M:%S"` COUNT LIBRARY SIZE. Processing $sample"
                countMap=$( cat $analysis_dir/mapping/hrm_developmental/$species/samtoolsFlagstat.$sample.remDup.txt | grep 'mapped (' | sed -E 's@\s\+\s.*@@g' )
                echo -e "$sample\t${countMap}" >> $analysis_dir/rpkm/hrm_developmental/mappedReadCounts.1.txt
            done
        done 1>>log/hrm_developmental/sortBam.picard_remDup.libSize.1.log

    #### 2. featureCounts
        script_featureCounts=/home/user/data2/lit/software/scripts/featureCounts.sh

        # for spe in human;do
        for spe in human rhesus mouse;do
            echo "***`date` FEATURECOUNTS. Processing $spe"
            SAMPLE=$(   
                ls $bamLst | 
                # head -n1 |
                while read bam;do
                sample=$( echo $bam | sed -E 's@^.*/@@g;s@.sorted.bam@@g' )
                printf "$analysis_dir/mapping/hrm_developmental/$spe/${sample}.remDup.bam|"
                done
            )

            gtf=$gtf_dir/${spe}.GRCh38.99.gtf

            bash $script_featureCounts \
            -g  $gtf  \
            -b  $SAMPLE \
            -o  $analysis_dir/featureCounts/hrm_developmental/$spe/featureCounts.${spe}.txt \
            -t  50 \
            -n  $spe \
            -ld log 1>log/hrm_developmental/featureCounts.$spe.log

    done
### ref ###
!


### 4. mapping qc ###

# 4.1 bam_stat before qc
. /home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh
conda activate RNA-seq
[ -d output/mapping_qc ] || mkdir -p output/mapping_qc
output_dir=output/mapping_qc

for i in sample;do for j in `seq 1 10`;do echo ${i}${j} ;done ;done |
# head -n 1 | \
while read sample;do
    #sorted bam has smaller size
    echo "***`date "+%Y-%m-%d %H:%M:%S"` RSeQC $sample.sorted.bam"
    time (bam_stat.py -i ${bam_dir}/$sample.bam > ${output_dir}/bam_stat_$sample.txt)
done


# 4.2 mark duplicates
script_remDup=/home/user/data2/lit/software/scripts/sambamba_remDup.sh

for i in sample;do for j in `seq 1 10`;do echo ${i}${j} ;done ;done |
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


for i in sample;do for j in `seq 1 10`;do echo ${i}${j} ;done ;done |
# head -n 1 | \
while read sample;do
    echo "***`date "+%Y-%m-%d %H:%M:%S"` Remove Duplicates and other qc for $sample.markDup.bam"
    time (samtools view -h -f 34 -F 1564 -@ 80  ${bam_dir}/${sample}.markDup.bam > ${bam_dir}/${sample}.clean.1.bam
    samtools view -h -f 18 -F 1580 -@ 80  ${bam_dir}/${sample}.markDup.bam > ${bam_dir}/${sample}.clean.2.bam
    samtools merge -@ 80 -f ${bam_dir}/${sample}.clean.bam ${bam_dir}/${sample}.clean.1.bam ${bam_dir}/${sample}.clean.2.bam 
    samtools flagstat -@ 80 ${bam_dir}/${sample}.clean.bam > ${output_dir}/samtoolsFlagstat.${sample}.clean.txt)
done 


# 4.4 bam_stat after qc

for i in sample;do for j in `seq 1 10`;do echo ${i}${j} ;done ;done |
# head -n 1 | \
while read sample;do
    #sorted bam has smaller size
    echo "***`date "+%Y-%m-%d %H:%M:%S"` RSeQC $sample.clean.bam"
    time (bam_stat.py -i ${bam_dir}/$sample.clean.bam > ${output_dir}/bam_stat_$sample.clean.txt)
done 


# nohup bash bin/run.sh >> log/mapping.qc.log 2>&1 &

### 4. mapping qc ###


### 5. featureCounts ###

[ -d output/featureCounts ] || mkdir output/featureCounts
output_dir=output/featureCounts
script_featureCounts=/home/user/data2/lit/software/scripts/featureCounts.sh

echo "*** featureCounts"

# seq 1 for test
SAMPLE=$( for i in sample;do for j in `seq 1 10`;do bam=/home/user/data/lit/project/6mA/RNA-seq/output/mapping/${i}${j}.clean.bam; 
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

# nohup bash bin/run.sh >> log/test.log 2>&1 &
# nohup bash bin/run.sh >> log/sort.mapping_qc.log 2>&1 &

<<'!'
### mapping quality control (samtools,sambamba) ###


    ### 1.remove low-quality or ambiguously mapped reads
    mapping_dir=/home/user/data/lit/project/6mA/mapping
    for i in 1 2;do
        # Read mapped in proper pair
        # remove Not primary alignment; Read fails platform/vendor quality checks; Supplementary alignment
        # remove Read unmapped; Mate unmapped
        samtools flagstat -@ 80 $mapping_dir/PC10_LY_MCC-M${i}.sam > qc/samtoolsFlagstat.PC10_LY_MCC-M${i}.unfil.txt
        samtools view -buSh -f 2 -F 2828 -q 20 -@ 80 $mapping_dir/PC10_LY_MCC-M${i}.sam > $mapping_dir/PC10_LY_MCC-M${i}.fil.bam
        samtools flagstat -@ 80 $mapping_dir/PC10_LY_MCC-M${i}.fil.bam > qc/samtoolsFlagstat.PC10_LY_MCC-M${i}.fil.txt
        samtools sort -O BAM -@ 80 $mapping_dir/PC10_LY_MCC-M${i}.fil.bam > $mapping_dir/PC10_LY_MCC-M${i}.sorted.bam
    done

    ### 2. remove duplicates

    for i in 1 2;do
        ## 2.1. Mark Duplicates
        echo -e "*** Mark Duplicates"
        time (sambamba markdup -t 10 $mapping_dir/PC10_LY_MCC-M${i}.sorted.bam $mapping_dir/PC10_LY_MCC-M${i}.markDup.bam)
        samtools flagstat -@ 80 $mapping_dir/PC10_LY_MCC-M${i}.markDup.bam > qc/samtoolsFlagstat.PC10_LY_MCC-M${i}.markDup.txt 

        ## 2.2 Remove Duplicates
        echo -e "*** Remove Duplicates"
        time (samtools view -b -F 1024 -@ 80 $mapping_dir/PC10_LY_MCC-M${i}.markDup.bam > $mapping_dir/PC10_LY_MCC-M${i}.remDup.bam)
        samtools flagstat -@ 80 $mapping_dir/PC10_LY_MCC-M${i}.remDup.bam > qc/samtoolsFlagstat.PC10_LY_MCC-M${i}.remDup.txt 
    done
### mapping quality control (samtools,sambamba) ###
!


multiqc -f  -o ./ output/mapping_qc/*clean.txt -n mapping_qc