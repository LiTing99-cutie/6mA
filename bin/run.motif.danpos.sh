################################################
#File Name: bin/run.motif.danpos.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Mon 14 Mar 2022 09:13:39 PM CST
################################################

#!/bin/sh 


set -e

# 0. get nucleosome or linker bed
    # nucleosome 50bp around nucleosome center; linker 50bp around linker center

    [ -d /home/user/data/lit/project/6mA/feature/iNPS/motif/danpos ] || mkdir -p /home/user/data/lit/project/6mA/feature/iNPS/motif/danpos
    [ -d log/iNPS/motif/danpos ] || mkdir -p log/iNPS/motif/danpos
    output_dir=/home/user/data/lit/project/6mA/feature/iNPS/motif/danpos

    genome_dir=/home/user/data/lit/database/public/genome/ce11/pengq/
    genome_fasta=${genome_dir}/Caenorhabditis_elegans.WBcel235.dna.fa
    genome_chrSize=$genome_dir/Caenorhabditis_elegans.WBcel235.dna.fa.size
    sort -k1,1 -k2,2n ${genome_chrSize} > $genome_dir/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size
    genome_sorted_chrSize=$genome_dir/Caenorhabditis_elegans.WBcel235.dna.fa.sorted.size
    ### nucleosome bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.nucleosome.bed"
        sed '1d' danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls | \
        awk -v OFS="\t" '{print $1,$2,$3,NR}' | sort -k1,1 -k2,2n > $output_dir/${sample}.nucleosome.bed 
        # only nucleosomes which have length >50 bp are kept
        less $output_dir/${sample}.nucleosome.bed | awk -v OFS="\t" '{if(($3-$2)>=50) print $0}' | \
        # only linkers which have length >50 bp are kept or the nucleosomes nearby are also removed
        sort -k1,1 -k2,2n | bedtools merge -d 50 -c 1,4 -o count,sum -i - | \
        awk -v OFS='\t' '$4==1{print $1,$2,$3,$5}' > $output_dir/${sample}.nucleosome.filter.bed
        # modify here
        # 50bp around nucleosome around center are kept for downstream analysis ( remove position which are negative and position which exceeds chromosome length )
        awk -v OFS="\t" '{print $1,(($2+$3)/2-25),(($2+$3)/2+25)}' $output_dir/${sample}.nucleosome.bed | grep -v '-' | \
        awk -v OFS='\t' '$1=="I" && $3<15072434{print $0} $1=="II" && $3<15279421{print $0} $1=="III" && $3<13783801{print $0} $1=="IV" && $3<17493829{print $0} $1=="V" && $3<20924180{print $0} $1=="X" && $3<17718942{print $0} $1=="MtDNA" && $3<13794{print $0}' > $output_dir/${sample}.nucleosome.center.bed
        nucleosome_number=`less $output_dir/${sample}.nucleosome.center.bed | wc -l`
        echo "${sample} nucleosome_number: ${nucleosome_number}"
    done )

    ### nucleosome bed ###
set +e
    ### linker bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.linker.bed"

        bedtools complement  -i $output_dir/${sample}.nucleosome.bed -g ${genome_sorted_chrSize} | awk -v OFS="\t" '{print $1,$2,$3,NR}' > \
        $output_dir/${sample}.linker.bed

        cut -f 4  $output_dir/${sample}.nucleosome.filter.bed | awk '{print $1+1}' > $output_dir/${sample}.nucleosome_id.add_1.txt
        cut -f 4  $output_dir/${sample}.nucleosome.filter.bed > $output_dir/${sample}.nucleosome_id.txt
        cat $output_dir/${sample}.nucleosome_id.add_1.txt $output_dir/${sample}.nucleosome_id.txt | sort -k1,1n | uniq > $output_dir/${sample}.linker_id.txt

        join -1 1 -2 4 $output_dir/${sample}.linker_id.txt  $output_dir/${sample}.linker.bed > $output_dir/${sample}.linker.filter.bed 
        
        awk -v OFS="\t" '{print $2,int(($3+$4)/2-25),int(($3+$4)/2+25)}' $output_dir/${sample}.linker.filter.bed | grep -v '-' > \
        $output_dir/${sample}.linker.center.bed
   
        linker_number=`less $output_dir/${sample}.linker.center.bed | wc -l`
        echo "${sample} linker_number: ${linker_number}"
    done )

    ### linker bed ###

set -e
# 1. get nucleosome or linker fasta

    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "get nucleosome fasta"
        # modify here 
        bedtools getfasta -fi $genome_fasta -bed $output_dir/${sample}.nucleosome.center.bed > \
                    $output_dir/${sample}.nucleosome.fasta &
        echo "get linker fasta"
        bedtools getfasta -fi $genome_fasta -bed $output_dir/${sample}.linker.center.bed > \
                    $output_dir/${sample}.linker.fasta &
    done 
    wait

# 2. extract motif
    # 3-base motif

    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    for j in nucleosome linker;do
        echo "extract motif for ${sample}.${j}"
        time ( cat $output_dir/${sample}.${j}.fasta | \
        # modify -p here
        seqkit locate -i -d -p NAN -j 20 | csvtk pretty -t | sed '1,2d' | awk -F '[\t:]' -v OFS='\t' '{print $1,$2}' | \
        awk -v OFS='\t' '{print $1,($6+$7)/2-1,($6+$7)/2,tolower($8)}' > ${output_dir}/${sample}.${j}.motif.bed )
    done
    done


# 3. stat motif feature

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
set +e
# nohup bash bin/run.motif.danpos.sh > log/iNPS/motif/danpos/part.0-3.log 2>&1 &