################################################
#File Name: run.danpos.n.1.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 25 Feb 2022 04:44:29 PM CST
################################################

#!/bin/sh 


################################################
#File Name: run.danpos.twoSamples.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 25 Feb 2022 01:20:58 PM CST
################################################

#!/bin/sh 


set -e

###  feature analysis ###

    [ -d /home/user/data/lit/project/6mA/feature/danpos/normed/method_1 ] || mkdir -p /home/user/data/lit/project/6mA/feature/danpos/normed/method_1
    [ -d log/danpos/normed/method_1 ] || mkdir -p log/danpos/normed/method_1
    output_dir=/home/user/data/lit/project/6mA/feature/danpos/normed/method_1
    log_dir=log/danpos/normed/method_1
    ce11_dir=/home/user/data/lit/database/public/genome/ce11/
    ce11_chrSize=$ce11_dir/ce11.fasta.size
    # awk -v OFS="\t" '{print $1,0,$2}' $ce11_chrSize > $ce11_dir/ce11.fasta.size.bed

    ### nucleosome bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.nucleosome.bed"
        less danpos/normed/${sample}/pooled/nucleosome_dyad.bed3 | \
        awk -v OFS="\t" '{print "chr"$1,$3-75,$3+75}' | grep -v '-' | \
        awk -v OFS='\t' '$1=="chrI" && $3<15072434{print $0} $1=="chrII" && $3<15279421{print $0} $1=="chrIII" && $3<13783801{print $0} $1=="chrIV" && $3<17493829{print $0} $1=="chrV" && $3<20924180{print $0} $1=="chrX" && $3<17718942{print $0} $1=="chrM" && $3<13794{print $0}' > $output_dir/${sample}.nucleosome.bed
    done )
    ### nucleosome bed ###

    ### linker bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.linker.bed"
        # alternative: bedtools complement
        bedtools subtract -a $ce11_dir/ce11.fasta.size.bed -b $output_dir/${sample}.nucleosome.bed > \
        $output_dir/${sample}.linker.0.bed

        i=0
        j=0
        cat $ce11_chrSize | \
        # head -n1 | \
        while read string;do 
            chr=$(echo -e ${string%[[:space:]]*})
            end=$(echo -e ${string#*[[:space:]]})
            let i=$i+1
            echo "processing $chr $i"
            sed "/$chr\t0/d;/$chr.*$end/d" $output_dir/${sample}.linker.$j.bed > \
            $output_dir/${sample}.linker.$i.bed
            let j=$j+1
        done

        cat $output_dir/${sample}.linker.7.bed > \
        $output_dir/${sample}.linker.bed 
        rm -rf $output_dir/${sample}.linker.[0-9].bed 
    done )


    ### linker bed ###


    ### 2. binned linker or nucleosome enrichment ###

    ### 2.1 bin ###
    time (
        for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "bin ${sample}"
        bedtools makewindows -b $output_dir/${sample}.nucleosome.bed  -n 10 -i winnum > \
        $output_dir/${sample}.nucleosome.binned.bed  2> $log_dir/${sample}.nucleosome.bin.log &
        bedtools makewindows -b $output_dir/${sample}.linker.bed  -n 10 -i winnum > \
        $output_dir/${sample}.linker.binned.bed 2> $log_dir/${sample}.linker.bin.log &

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
        ' $output_dir/${sample}.nucleosome.binned.bed > $output_dir/${sample}.nucleosome.binned.reorg.bed

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
        ' $output_dir/${sample}.linker.binned.bed > $output_dir/${sample}.linker.binned.reorg.bed 
        done )

    ### 2.1 bin ###

    ### 2.2 6mA level ###
        [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt
        time (
        for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
            for j in nucleosome linker;do
                for i in 1 2 3 4 5;do
                    echo "*** process ${sample}.$j.bin_$i"
                    less $output_dir/${sample}.$j.binned.reorg.bed | \
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
                for i in 1 2 3 4 5;do
                    mA=`less $output_dir/${sample}.$j.bin_$i.intersect.bed | wc -l`
                    A=`grep -i -o a $output_dir/${sample}.$j.bin_$i.fasta | wc -l`
                    FreqSum=`less $output_dir/${sample}.$j.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`
                    echo -e "${sample}.$j.bin_$i\t${mA}\t${A}\t${FreqSum}" >> $output_dir/binned.level.txt
                done
            done
        done
    ### 2.2 6mA level ###

    ### 2. binned linker or nucleosome enrichment ###

### feature analysis ###

# 2022-02-25
# [ -d log/danpos/normed/method_1 ] || mkdir -p log/danpos/normed/method_1
# nohup bash run.danpos.n.1.sh > log/danpos/normed/method_1/run.log 2>&1 &