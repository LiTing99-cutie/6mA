################################################
#File Name: run.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 08 Feb 2022 08:47:38 PM CST
################################################

#!/bin/sh 

set -e

<<'!'
### 2022-02-08 build bowtie indexes (bowtie2) ###
    [ -d log ] || mkdir log 
    [ -d index ] || mkdir index 
    [ -d /home/user/data/lit/project/6mA/mapping ] || mkdir /home/user/data/lit/project/6mA/mapping
    [ -d qc ] || mkdir qc
    genome_fasta=/home/user/data/lit/database/public/genome/ce11/ce11.fa.gz

    pushd index
    time (bowtie2-build --threads 20 $genome_fasta ce11) 
    popd 

    # nohup bash run.sh >log/bowtie2-build.log 2>&1 &
!

<<'!'
### 2022-02-09 rawdata quality control (fastp) ###
    data_dir=/home/user/data/lit/project/6mA/data/MNase/
    . "/home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh"
    conda activate RNA-seq
    for i in 1 2;do
        fastp --thread=8 \
        -i $data_dir/PC10_LY_MCC-M${i}/PC10_LY_MCC-M${i}.1.clean.fq.gz -o $data_dir/PC10_LY_MCC-M${i}.1.fastq.gz \
        -I $data_dir/PC10_LY_MCC-M${i}/PC10_LY_MCC-M${i}.2.clean.fq.gz -O $data_dir/PC10_LY_MCC-M${i}.2.fastq.gz \
        -f 5 \
        -F 5 \
        -h qc/PC10_LY_MCC-M${i}.html \
        -j qc/PC10_LY_MCC-M${i}.json
    done

    # nohup bash run.sh >log/fastp.log 2>&1 &
!

<<'!'
### 2022-02-09 mapping (bowtie2) ###
    data_dir=/home/user/data/lit/project/6mA/data/MNase/
    mapping_dir=/home/user/data/lit/project/6mA/mapping
    for i in 1 2;do
        fq_1=$data_dir/PC10_LY_MCC-M${i}.1.fastq.gz
        fq_2=$data_dir/PC10_LY_MCC-M${i}.2.fastq.gz
        time (bowtie2 -p 20 -x index/ce11 -1 $fq_1 -2 $fq_2 -S $mapping_dir/PC10_LY_MCC-M${i}.sam)
    done

    # nohup bash run.sh >log/bowtie2.mapping.log 2>&1 &
!


<<'!'
### 2022-02-09 mapping quality control (samtools,sambamba) ###

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

    # nohup bash run.sh >log/mapping.qc.log 2>&1 &
!

<<'!'
### 2022-02-10 bam to bed (bedtools) ###
    mapping_dir=/home/user/data/lit/project/6mA/mapping

    for i in 1 2;do
        bedtools bamtobed -i $mapping_dir/PC10_LY_MCC-M${i}.remDup.bam > $mapping_dir/PC10_LY_MCC-M${i}.bed
    done

    # nohup bash run.sh >log/bamtoBed.log 2>&1 &
!

<<'!'
### 2022-02-10 iNPS ###
    # output_dir=/home/user/data/lit/project/6mA/iNPS_output
    # mapping_dir=/home/user/data/lit/project/6mA/mapping
    # [ -d /home/user/data/lit/project/6mA/iNPS_output ] || mkdir -p /home/user/data/lit/project/6mA/iNPS_output
     
    # for i in 1 2;do
    #     samtools sort -n -O BAM -@ 80 $mapping_dir/PC10_LY_MCC-M${i}.remDup.bam > $mapping_dir/PC10_LY_MCC-M${i}.remDup.sorted.bam
    #     bedtools bamtobed -bedpe -i $mapping_dir/PC10_LY_MCC-M${i}.remDup.sorted.bam > $mapping_dir/PC10_LY_MCC-M${i}.bedpe
    #     awk  -v OFS='\t' '{print $1,$2,$6}'  $mapping_dir/PC10_LY_MCC-M${i}.bedpe > $mapping_dir/PC10_LY_MCC-M${i}.bedpe.bed
    #     iNPS=/home/user/data2/lit/software/bin/iNPS_V1.2.2.py
    #     time (python3 $iNPS --s_p=p -i $mapping_dir/PC10_LY_MCC-M${i}.bedpe.bed -o $output_dir/PC10_LY_MCC-M${i})
    # done

    # nohup bash run.sh >log/iNPS.log 2>&1 &
    

    ### wrong ###
        # output_dir=/home/user/data/lit/project/6mA/iNPS_output
        # mapping_dir=/home/user/data/lit/project/6mA/mapping
        # [ -d /home/user/data/lit/project/6mA/iNPS_output ] || mkdir -p /home/user/data/lit/project/6mA/iNPS_output
        # for i in 1 2;do
        #     awk  -v OFS='\t' '{print $1,$2,$3}'  $mapping_dir/PC10_LY_MCC-M${i}.bed > $mapping_dir/PC10_LY_MCC-M${i}.input.bed
        #     iNPS=/home/user/data2/lit/software/bin/iNPS_V1.2.2.py
        #     time (python3 $iNPS --s_p=p -i $mapping_dir/PC10_LY_MCC-M${i}.input.bed -o $output_dir/PC10_LY_MCC-M${i}.wrong)
        # done
        # nohup bash run.sh >log/iNPS.wrong.log 2>&1 &
    ### wrong ###

    ### test ###
        # iNPS=/home/user/data2/lit/software/bin/iNPS_V1.2.2.py
        # output_dir=/home/user/data/lit/project/6mA/iNPS_output
        # test_bed=/home/user/data/lit/project/6mA/GSM849959_GA2807_CMT1_shH2A.Z-2d_MNase_0.1U_r520l2.bed
        # time (python3 $iNPS --s_p=p -i ${test_bed} -o $output_dir/test)

        # # nohup bash run.sh >log/iNPS.test.log 2>&1 &
    ### test ###

### 2022-02-10 iNPS ###
!

<<'!'
### 2022-02-11 visualization ###



    output_dir=/home/user/data/lit/project/6mA/iNPS_output

    ce11_dir=/home/user/data/lit/database/public/genome/ce11/
    # gunzip -c $ce11_dir/ce11.fa.gz > $ce11_dir/ce11.fa
    # samtools faidx $ce11_dir/ce11.fa && \
    # awk '{print $1"\t"$2}' $ce11_dir/ce11.fa.fai > $ce11_dir/ce11.fasta.size 
    ce11_chrSize=$ce11_dir/ce11.fasta.size 

    # awk -v OFS="\t" '{print $1,$2,$3,$6}' $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_bed | tail -n +3 | \
    # sort -k1,1 -k2,2n \
    # > $output_dir/PC10_LY_MCC-M1.bedpe_chrI.bdg

    # bedGraphToBigWig $output_dir/PC10_LY_MCC-M1.bedpe_chrI.bdg $ce11_chrSize $output_dir/PC10_LY_MCC-M1.bedpe_chrI.peak.bigwig

    # for i in 1 2;do
    #     for chr in chrI chrII chrIII chrIV chrV chrX chrM; do

    awk -v OFS="\t" 'NR==2 {sub(/chromosome/,"chrom"); print $0}' PC10_LY_MCC-M1.bedpe_chrI.like_wig | sed -n '2p' > \
    $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig.title

    awk -v OFS="\t" '{print $1,$2}' $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig | sed '1,3d' > \
    $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig.titleRm.ori

    cat $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig.title $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig.titleRm.ori > \
    $output_dir/PC10_LY_MCC-M1.bedpe_chrI.OriginalProfile.wig
    
    # wigToBigWig $output_dir/PC10_LY_MCC-M1.bedpe_chrI.OriginalProfile.wig $ce11_chrSize \
    # $output_dir/PC10_LY_MCC-M1.bedpe_chrI.OriginalProfile.bigwig

    awk -v OFS="\t" '{print $1,$3}' $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig | sed '1,3d' > \
    $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig.titleRm.smoo

    cat $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig.title $output_dir/PC10_LY_MCC-M1.bedpe_chrI.like_wig.titleRm.smoo > \
    $output_dir/PC10_LY_MCC-M1.bedpe_chrI.SmoothProfile.wig

    # wigToBigWig $output_dir/PC10_LY_MCC-M1.bedpe_chrI.OriginalProfile.wig $ce11_chrSize \
    # $output_dir/PC10_LY_MCC-M1.bedpe_chrI.SmoothProfile.bigwig

    # sudo ln -s $output_dir/PC10_LY_MCC-M1.bedpe_chrI.peak.bigwig \
    # /home/user/data2/rbase/tomcat_base/webapps/RhesusBase/lit/PC10_LY_MCC-M1.bedpe_chrI.peak.bigwig

    # sudo ln -s $output_dir/PC10_LY_MCC-M1.bedpe_chrI.SmoothProfile.bigwig \
    # /home/user/data2/rbase/tomcat_base/webapps/RhesusBase/lit/PC10_LY_MCC-M1.bedpe_chrI.SmoothProfile.bigwig

    # sudo ln -s $output_dir/PC10_LY_MCC-M1.bedpe_chrI.OriginalProfile.bigwig \
    # /home/user/data2/rbase/tomcat_base/webapps/RhesusBase/lit/PC10_LY_MCC-M1.bedpe_chrI.OriginalProfile.bigwig

    # nohup bash run.sh >log/visualization.log 2>&1 &

### 2022-02-14 visualization ###
!


### 2022-02-16 feature analysis ###


    output_dir=/home/user/data/lit/project/6mA/iNPS_output
    ce11_dir=/home/user/data/lit/database/public/genome/ce11/
    ce11_chrSize=$ce11_dir/ce11.fasta.size
    # awk -v OFS="\t" '{print $1,0,$2}' $ce11_chrSize > $ce11_dir/ce11.fasta.size.bed

    ### nucleosome bed ###
    time (
    for i in 1 2;do
        echo "*** generating PC10_LY_MCC-M${i}_Gathering.nucleosome.bed"
        sed '1,8d' $output_dir/PC10_LY_MCC-M${i}_Gathering.like_bed | \
        awk -v OFS="\t" '{print $1,(($2+$3+10)/2-73),(($2+$3+10)/2+73)}' | grep -v '-' > \
        PC10_LY_MCC-M${i}_Gathering.nucleosome.bed
    done )
    ### nucleosome bed ###


    ### linker bed ###
    time (
    for id in 1 2;do
        echo "*** generating PC10_LY_MCC-M${id}_Gathering.linker.bed"
        # alternative: bedtools complement
        bedtools subtract -a $ce11_dir/ce11.fasta.size.bed -b PC10_LY_MCC-M${id}_Gathering.nucleosome.bed > \
        $output_dir/PC10_LY_MCC-M${id}_Gathering.linker.0.bed

        i=0
        j=0
        cat $ce11_chrSize | \
        # head -n1 | \
        while read string;do 
            chr=$(echo -e ${string%[[:space:]]*})
            end=$(echo -e ${string#*[[:space:]]})
            let i=$i+1
            echo "processing $chr $i"
            sed "/$chr\t0/d;/$chr.*$end/d" $output_dir/PC10_LY_MCC-M${id}_Gathering.linker.$j.bed > \
            $output_dir/PC10_LY_MCC-M${id}_Gathering.linker.$i.bed
            let j=$j+1
        done

        cat $output_dir/PC10_LY_MCC-M${id}_Gathering.linker.7.bed > \
        PC10_LY_MCC-M${id}_Gathering.linker.bed 
        rm -rf $output_dir/PC10_LY_MCC-M${id}_Gathering.linker.[0-9].bed 
    done )

    ### linker bed ###


    ### 1. linker or nucleosome enrichment ###
    [ -f level.txt ] && rm -rf level.txt
    time (
    for i in 1 2; do
        for j in nucleosome linker; do
            echo "processing PC10_LY_MCC-M${i}.${j}"
            bedtools intersect -wb -a PC10_LY_MCC-M${i}_Gathering.${j}.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
            awk -v OFS="\t" '{print $4,$5,$6,$7}' > PC10_LY_MCC-M${i}_Gathering.${j}.intersect.bed 

            bedtools getfasta -fi $ce11_dir/ce11.fa -bed PC10_LY_MCC-M${i}_Gathering.${j}.bed > \
            $output_dir/PC10_LY_MCC-M${i}_Gathering.${j}.fasta 

            FreqSum=`less PC10_LY_MCC-M${i}_Gathering.${j}.intersect.bed | awk '{sum += $4} END {print sum}'`
            mA=`less PC10_LY_MCC-M${i}_Gathering.${j}.intersect.bed | wc -l`
            A=`grep -i -o a $output_dir/PC10_LY_MCC-M${i}_Gathering.${j}.fasta | wc -l`
            # grep -v '^>' $output_dir/PC10_LY_MCC-M${i}_Gathering.${j}.fasta | grep -i -o [a-z] | wc -l
            echo -e "PC10_LY_MCC-M${i}.${j}\t${mA}\t${A}\t${FreqSum}" >> level.txt
        done 
    done )

    ### 1. linker or nucleosome enrichment ###



    ### 2. binned linker or nucleosome enrichment ###

    ### 2.1 bin ###
    time (
    for i in 1 2;do
        echo "bin PC10_LY_MCC-M${i}"
        bedtools makewindows -b PC10_LY_MCC-M${i}_Gathering.nucleosome.bed  -n 10 -i winnum > \
        $output_dir/PC10_LY_MCC-M${i}_Gathering.nucleosome.binned.bed  2> log/PC10_LY_MCC-M${i}.nucleosome.bin.log 
        bedtools makewindows -b PC10_LY_MCC-M${i}_Gathering.linker.bed  -n 10 -i winnum > \
        $output_dir/PC10_LY_MCC-M${i}_Gathering.linker.binned.bed 2>log/PC10_LY_MCC-M${i}.linker.bin.log 

        awk '$4==1 {print $1,$2,$3,$4+4}
        $4==2 {print $1,$2,$3,$4+2}
        $4==3 {print $1,$2,$3,$4}
        $4==4 {print $1,$2,$3,$4-2}
        $4==5 {print $1,$2,$3,$4-4}
        $4==6 {print $1,$2,$3,$4-5}
        $4==7 {print $1,$2,$3,$4-5}
        $4==8 {print $1,$2,$3,$4-5}
        $4==9 {print $1,$2,$3,$4-5}
        $4==10 {print $1,$2,$3,$4-5}
        ' $output_dir/PC10_LY_MCC-M${i}_Gathering.nucleosome.binned.bed > $output_dir/PC10_LY_MCC-M${i}_Gathering.nucleosome.binned.reorg.bed

        awk '$4==1 {print $1,$2,$3,$4}
        $4==2 {print $1,$2,$3,$4}
        $4==3 {print $1,$2,$3,$4}
        $4==4 {print $1,$2,$3,$4}
        $4==5 {print $1,$2,$3,$4}
        $4==6 {print $1,$2,$3,$4-1}
        $4==7 {print $1,$2,$3,$4-3}
        $4==8 {print $1,$2,$3,$4-5}
        $4==9 {print $1,$2,$3,$4-7}
        $4==10 {print $1,$2,$3,$4-9}
        ' $output_dir/PC10_LY_MCC-M${i}_Gathering.linker.binned.bed > $output_dir/PC10_LY_MCC-M${i}_Gathering.linker.binned.reorg.bed
    done )

    ### 2.1 bin ###

    ### 2.2 6mA level ###
        [ -f binned.level.txt ] && rm -rf binned.level.txt
        time (
        # for i in 1;do
        for i in 1 2 3 4 5;do
            for j in nucleosome linker;do
            # for j in nucleosome;do
                for id in 1 2;do
                # for id in 1;do
                    echo "*** process PC10_LY_MCC-M${id}.$j.bin_$i"
                    less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.binned.reorg.bed | \
                            grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
                                $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed
                    bedtools intersect -wb -a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                    awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.bed 
                    bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.bed > \
                    $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.fasta
                    mA=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.bed | wc -l`
                    A=`grep -i -o a $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.fasta | wc -l`
                    FreqSum=`less $output_dir/PC10_LY_MCC-M${id}_Gathering.$j.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`
                    echo -e "PC10_LY_MCC-M${id}.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> binned.level.txt
                done
            done
        done ) 
    ### 2.2 6mA level ###

    ### 2. binned linker or nucleosome enrichment ###

### 2022-02-16 feature analysis ###

# nohup bash run.sh >log/feature.1.log 2>&1 &


