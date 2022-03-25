################################################
#File Name: run.danpos.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 24 Feb 2022 01:27:24 PM CST
################################################

#!/bin/sh 

set -e

### 2022-02-24 feature analysis ###

    [ -d /home/user/data/lit/project/6mA/feature/danpos/method_1 ] || mkdir -p /home/user/data/lit/project/6mA/feature/danpos/method_1
    [ -d log/danpos ] || mkdir -p log/danpos
    output_dir=/home/user/data/lit/project/6mA/feature/danpos/method_1
    log_dir=log/danpos
    ce11_dir=/home/user/data/lit/database/public/genome/ce11/
    ce11_chrSize=$ce11_dir/ce11.fasta.size
    # awk -v OFS="\t" '{print $1,0,$2}' $ce11_chrSize > $ce11_dir/ce11.fasta.size.bed

    ### nucleosome bed ###
    time (
        echo "*** generating nucleosome.bed"
        less pooled/nucleosome_dyad.bed3 | \
        awk -v OFS="\t" '{print "chr"$1,$3-73,$3+73}' | grep -v '-' | \
        awk -v OFS='\t' '$1=="chrI" && $3<15072434{print $0} $1=="chrII" && $3<15279421{print $0} $1=="chrIII" && $3<13783801{print $0} $1=="chrIV" && $3<17493829{print $0} $1=="chrV" && $3<20924180{print $0} $1=="chrX" && $3<17718942{print $0} $1=="chrM" && $3<13794{print $0}' > $output_dir/danpos.nucleosome.bed)
    ### nucleosome bed ###

    ### linker bed ###
    time (
        echo "*** generating linker.bed"
        # alternative: bedtools complement
        bedtools subtract -a $ce11_dir/ce11.fasta.size.bed -b $output_dir/danpos.nucleosome.bed > \
        $output_dir/danpos.linker.0.bed

        i=0
        j=0
        cat $ce11_chrSize | \
        # head -n1 | \
        while read string;do 
            chr=$(echo -e ${string%[[:space:]]*})
            end=$(echo -e ${string#*[[:space:]]})
            let i=$i+1
            echo "processing $chr $i"
            sed "/$chr\t0/d;/$chr.*$end/d" $output_dir/danpos.linker.$j.bed > \
            $output_dir/danpos.linker.$i.bed
            let j=$j+1
        done

        cat $output_dir/danpos.linker.7.bed > \
        $output_dir/danpos.linker.bed 
        rm -rf $output_dir/danpos.linker.[0-9].bed )

    ### linker bed ###


    ### 1. linker or nucleosome enrichment ###
    [ -f $output_dir/level.txt ] && rm -rf $output_dir/level.txt
    time (
    for j in nucleosome linker; do
            echo "processing danpos.$j.bed "
            bedtools intersect -wb -a $output_dir/danpos.${j}.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
            awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/danpos.${j}.intersect.bed 

            bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/danpos.${j}.bed > \
            $output_dir/danpos.${j}.fasta 

            FreqSum=`less $output_dir/danpos.${j}.intersect.bed | awk '{sum += $4} END {print sum}'`
            mA=`less $output_dir/danpos.${j}.intersect.bed | wc -l`
            A=`grep -i -o a $output_dir/danpos.${j}.fasta | wc -l`
            echo -e "danpos.${j}\t${mA}\t${A}\t${FreqSum}" >> $output_dir/level.txt
    done )

    ### 1. linker or nucleosome enrichment ###


    ### 2. binned linker or nucleosome enrichment ###

    ### 2.1 bin ###
    time (
        echo "bin danpos"
        bedtools makewindows -b $output_dir/danpos.nucleosome.bed  -n 10 -i winnum > \
        $output_dir/danpos.nucleosome.binned.bed  2> $log_dir/danpos.nucleosome.bin.log &
        bedtools makewindows -b $output_dir/danpos.linker.bed  -n 10 -i winnum > \
        $output_dir/danpos.linker.binned.bed 2> $log_dir/danpos.linker.bin.log &

        wait

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
        ' $output_dir/danpos.nucleosome.binned.bed > $output_dir/danpos.nucleosome.binned.reorg.bed

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
        ' $output_dir/danpos.linker.binned.bed > $output_dir/danpos.linker.binned.reorg.bed )

    ### 2.1 bin ###

    ### 2.2 6mA level ###
        [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt
        time (
        for j in nucleosome linker;do
            for i in 1 2 3 4 5;do
                    echo "*** process danpos.$j.bin_$i"
                    less $output_dir/danpos.$j.binned.reorg.bed | \
                            grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > \
                                $output_dir/danpos.$j.bin_$i.bed
                    bedtools intersect -wb -a $output_dir/danpos.$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                    awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/danpos.$j.bin_$i.intersect.bed &

                    bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/danpos.$j.bin_$i.bed > \
                    $output_dir/danpos.$j.bin_$i.fasta &
            done
            wait
        done ) 


        for i in 1 2 3 4 5;do
            for j in nucleosome linker;do
                    mA=`less $output_dir/danpos.$j.bin_$i.intersect.bed | wc -l`
                    A=`grep -i -o a $output_dir/danpos.$j.bin_$i.fasta | wc -l`
                    FreqSum=`less $output_dir/danpos.$j.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`
                    echo -e "danpos.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> $output_dir/binned.level.txt
            done
        done
    ### 2.2 6mA level ###

    ### 2. binned linker or nucleosome enrichment ###

### 2022-02-24 feature analysis ###

# 2022-02-24

# nohup bash run.danpos.sh > log/danpos/run.log 2>&1 &
