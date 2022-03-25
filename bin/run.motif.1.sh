################################################
#File Name: bin/run.motif.1.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 14 Mar 2022 11:34:59 AM CST
################################################

#!/bin/sh 

set -e
<<'!'
# 0. get nucleosome or linker bed
# nucleosome 50bp around nucleosome center; linker 50bp around linker center
    output_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif
    ce11_dir=/home/user/data/lit/database/public/genome/ce11/
    ce11_chrSize=$ce11_dir/ce11.fasta.size
    ### nucleosome bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.nucleosome.bed"
        sed '1,8d' /home/user/data/lit/project/6mA/iNPS_output/${sample}_Gathering.like_bed | \
        # only nucleosomes which have length >50 bp are kept
        awk -v OFS="\t" '{if($5>=50) print $1,$2,$3}' | \
        # only linkers which have length >50 bp are kept or the nucleosomes nearby are also removed
        sort -k1,1 -k2,2n | bedtools merge -d 50 -c 1 -o count -i - | \
        awk -v OFS='\t' '$4==1{print $1,$2,$3}' > $output_dir/${sample}.nucleosome.bed
        # modify here
        # 50bp around nucleosome around center are kept for downstream analysis ( remove position which are negative and position which exceeds chromosome length )
        awk -v OFS="\t" '{print $1,(($2+$3+10)/2-25),(($2+$3+10)/2+25)}' $output_dir/${sample}.nucleosome.bed | grep -v '-' | \
        awk -v OFS='\t' '$1=="chrI" && $3<15072434{print $0} $1=="chrII" && $3<15279421{print $0} $1=="chrIII" && $3<13783801{print $0} $1=="chrIV" && $3<17493829{print $0} $1=="chrV" && $3<20924180{print $0} $1=="chrX" && $3<17718942{print $0} $1=="chrM" && $3<13794{print $0}' > $output_dir/${sample}.nucleosome.center.bed
        nucleosome_number=`less $output_dir/${sample}.nucleosome.bed | wc -l`
        echo "${sample} nucleosome_number: ${nucleosome_number}"
    done )

    ### nucleosome bed ###

    ### linker bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.linker.bed"
        # alternative: bedtools complement
        bedtools subtract -a $ce11_dir/ce11.fasta.size.bed -b $output_dir/${sample}.nucleosome.bed > \
        $output_dir/${sample}.linker.0.bed

        # remove the first linker and the last linker
        i=0
        j=0
        while read string;do 
            chr=$(echo -e ${string%[[:space:]]*})
            end=$(echo -e ${string#*[[:space:]]})
            let i=$i+1
            echo "processing $chr $i"
            sed "/$chr\t0/d;/$chr.*$end/d" $output_dir/${sample}.linker.$j.bed > \
            $output_dir/${sample}.linker.$i.bed
            let j=$j+1
        done < $ce11_chrSize

        echo "*** processing ${sample}.linker.${i}.bed"
        cat $output_dir/${sample}.linker.$i.bed | awk -v OFS="\t" '{print $1,int(($2+$3)/2-25),int(($2+$3)/2+25)}' > \
        $output_dir/${sample}.linker.bed 
        rm -rf $output_dir/${sample}.linker.[0-9].bed 

        linker_number=`less $output_dir/${sample}.linker.bed | wc -l`
        echo "${sample} linker_number: ${linker_number}"
    done )

    ### linker bed ###

    # nohup bash bin/run.motif.1.sh > log/iNPS/motif/get_bed.log 2>&1 &
!

<<'!'
# 1. get nucleosome or linker fasta
    [ -d /home/user/data/lit/project/6mA/feature/iNPS/motif ] || mkdir -p /home/user/data/lit/project/6mA/feature/iNPS/motif
    [ -d log/iNPS/motif ] || mkdir -p log/iNPS/motif
    input_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif
    output_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif
    ce11_dir=/home/user/data/lit/database/public/genome/ce11/

    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "get nucleosome fasta"
        bedtools getfasta -fi $ce11_dir/ce11.fa -bed $input_dir/${sample}.nucleosome.bed > \
                    $input_dir/${sample}.nucleosome.fasta &
        echo "get linker fasta"
        bedtools getfasta -fi $ce11_dir/ce11.fa -bed $input_dir/${sample}.linker.bed > \
                    $input_dir/${sample}.linker.fasta &
    done 
    wait

# 2. extract motif
# 3-base motif
    input_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif
    output_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif

    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    for j in nucleosome linker;do
        echo "extract motif for ${sample}.${j}"
        time ( cat $input_dir/${sample}.${j}.fasta | \
        # modify -p here
        seqkit locate -i -d -p NAN -j 20 | csvtk pretty -t | sed '1,2d' | awk -F '[\t:]' -v OFS='\t' '{print $1,$2}' | \
        awk -v OFS='\t' '{print $1,($6+$7)/2-1,($6+$7)/2,tolower($8)}' > ${output_dir}/${sample}.${j}.motif.bed )
    done
    done

# nohup bash bin/run.motif.1.sh > log/iNPS/motif/part.1-2.log 2>&1 &
!

# 3. stat motif feature

    
    input_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif
    output_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif
    [ -f $output_dir/motif.stat.txt ] && rm -rf $output_dir/motif.stat.txt

    time ( for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    for j in nucleosome linker;do
        echo "stat for ${sample}.${j}"
        awk '{print $4}' $output_dir/${sample}.${j}.motif.bed | sort | uniq -i > $output_dir/${sample}.${j}.motif.txt
        cat $output_dir/${sample}.${j}.motif.txt | \
        while read motif;do 
            echo "processing $motif"
            awk '{if($4=="'$motif'") print$0}' $output_dir/${sample}.${j}.motif.bed > $output_dir/${sample}.${j}.${motif}.bed
            motif_n=`less ${output_dir}/${sample}.${j}.${motif}.bed | wc -l`
            echo -e "${sample}\t${j}\t${motif}\t${motif_n}" >> $output_dir/motif.stat.txt 
        done
    done
    done )

# nohup bash bin/run.motif.1.sh > log/iNPS/motif/part.3.log 2>&1 &
