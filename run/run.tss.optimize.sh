################################################
#File Name: run.tss.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 23 Feb 2022 05:08:03 PM CST
################################################

#!/bin/sh 




### 3. binned TSS enrichment ###


set -e

# output_dir=/home/user/data/lit/project/6mA/iNPS_output/
# ce11_dir=/home/user/data/lit/database/public/genome/ce11/
# ce11_chrSize=$ce11_dir/ce11.fasta.size
# gene_dir=/home/user/data/lit/database/public/annotation/gene_and_gene_predictions/


### 3.1 bin ###

    ### bin tss 2kb upstream and downstream ###
        # grep -v '^#' $gene_dir/ce11.refGene.bed | awk -v OFS='\t' '{print $3,$5-2000,$5+2000}' |  grep -v '-' | \
        # awk -v OFS='\t' '$1=="chrI" && $3<15072434{print $0} $1=="chrII" && $3<15279421{print $0} $1=="chrIII" && $3<13783801{print $0} $1=="chrIV" && $3<17493829{print $0} $1=="chrV" && $3<20924180{print $0} $1=="chrX" && $3<17718942{print $0} $1=="chrM" && $3<13794{print $0}' > $gene_dir/ce11.tss.bed
        
        # nohup bash -c "time (less $gene_dir/ce11.tss.bed | bedtools makewindows -b - -w 20 -i winnum > $output_dir/ce11.tss.binned.20.bed)" &
    ### bin tss 2kb upstream and downstream ###

### 3.1 bin ###


### like_wiggle to wig ###
    # done in run.tss.sh
### like_wiggle to wig ###


### wig to bed ###

    # time ( 
    # for i in 1 2;do
    #     for j in chrI chrII chrIII chrIV chrV chrX chrM; do
    #         for profile in OriginalProfile SmoothProfile;do
    #             echo "wig2bed PC10_LY_MCC-M${i}_$j.$profile"
    #             wig2bed < $output_dir/PC10_LY_MCC-M${i}_$j.$profile.wig | awk -v OFS="\t" '{print $1,$2,$3,$5}' > \
    #             $output_dir/PC10_LY_MCC-M${i}_$j.$profile.bed
    #         done
    #     done
    # done ) 

    # done in run.tss.sh

### wig to bed ###

### nucleosome wigbed: mutiple chromosome to one chromosome ###
    # time ( for i in 1 2;do
    #     for profile in OriginalProfile SmoothProfile;do
    #         for j in chrI chrII chrIII chrIV chrV chrX chrM; do
    #             cat $output_dir/PC10_LY_MCC-M${i}_$j.$profile.bed >> $output_dir/PC10_LY_MCC-M${i}.$profile.combined.bed
    #         done
    #     done
    # done ) &
### nucleosome wigbed: mutiple chromosome to one chromosome ###


### 2022-02-24 change output dir and manual multi-threading ###

### change output dir ###

[ -d /home/user/data/lit/project/6mA/feature/tss ] || mkdir -p /home/user/data/lit/project/6mA/feature/tss

[ -d log/tss ] || mkdir -p log/tss

output_dir=/home/user/data/lit/project/6mA/feature/tss
ce11_dir=/home/user/data/lit/database/public/genome/ce11/
ce11_chrSize=$ce11_dir/ce11.fasta.size
gene_dir=/home/user/data/lit/database/public/annotation/gene_and_gene_predictions/

###  3.2 6mA level and nucleosome occupancy ###
    [ -f $output_dir/tss.mAlevel.occu.txt ] && rm -rf $output_dir/tss.mAlevel.occu.txt

    ### manual multi-threading ###

    ### block 1 ###
        time ( 
        for j in ce11.tss;do
            for i in `seq 1 50`;do
                echo "*** process $j.bin_$i"
                less /home/user/data/lit/project/6mA/iNPS_output/$j.binned.20.bed | \
                # head -n 40000 | \
                grep -P "\t$i$" | awk -v OFS='\t' '{print $1,$2,$3}' > $output_dir/$j.bin_$i.bed

                bedtools intersect -wb -a $output_dir/$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/$j.bin_$i.intersect.bed &

                bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/$j.bin_$i.bed > $output_dir/$j.bin_$i.fasta &

                bedtools intersect -wo -a /home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M1.SmoothProfile.combined.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M1.bed &
                
                bedtools intersect -wo -a /home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M2.SmoothProfile.combined.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M2.bed &

            done 
            wait
        done ) 
    ### block 1 ###

    ### block 2 ###
        time ( 
        for j in ce11.tss;do
            for i in `seq 51 100`;do
                echo "*** process $j.bin_$i"
                less /home/user/data/lit/project/6mA/iNPS_output/$j.binned.20.bed | \
                # head -n 40000 | \
                grep -P "\t$i$" | awk -v OFS='\t' '{print $1,$2,$3}' > $output_dir/$j.bin_$i.bed

                bedtools intersect -wb -a $output_dir/$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/$j.bin_$i.intersect.bed &

                bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/$j.bin_$i.bed > $output_dir/$j.bin_$i.fasta &

                bedtools intersect -wo -a /home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M1.SmoothProfile.combined.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M1.bed &
                
                bedtools intersect -wo -a /home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M2.SmoothProfile.combined.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M2.bed &

            done 
            wait
        done ) 
    ### block 2 ###   

    ### block 3 ###
        time ( 
        for j in ce11.tss;do
            for i in `seq 101 150`;do
                echo "*** process $j.bin_$i"
                less /home/user/data/lit/project/6mA/iNPS_output/$j.binned.20.bed | \
                # head -n 40000 | \
                grep -P "\t$i$" | awk -v OFS='\t' '{print $1,$2,$3}' > $output_dir/$j.bin_$i.bed

                bedtools intersect -wb -a $output_dir/$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/$j.bin_$i.intersect.bed &

                bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/$j.bin_$i.bed > $output_dir/$j.bin_$i.fasta &

                bedtools intersect -wo -a /home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M1.SmoothProfile.combined.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M1.bed &
                
                bedtools intersect -wo -a /home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M2.SmoothProfile.combined.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M2.bed &

            done 
            wait
        done ) 
    ### block 3 ###  

    ### block 4 ###
        time ( 
        for j in ce11.tss;do
            for i in `seq 151 200`;do
                echo "*** process $j.bin_$i"
                less /home/user/data/lit/project/6mA/iNPS_output/$j.binned.20.bed | \
                # head -n 40000 | \
                grep -P "\t$i$" | awk -v OFS='\t' '{print $1,$2,$3}' > $output_dir/$j.bin_$i.bed

                bedtools intersect -wb -a $output_dir/$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
                awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/$j.bin_$i.intersect.bed &

                bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/$j.bin_$i.bed > $output_dir/$j.bin_$i.fasta &

                bedtools intersect -wo -a /home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M1.SmoothProfile.combined.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M1.bed &
                
                bedtools intersect -wo -a /home/user/data/lit/project/6mA/iNPS_output/PC10_LY_MCC-M2.SmoothProfile.combined.bed -b $output_dir/$j.bin_$i.bed \
                > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M2.bed &

            done 
            wait
        done ) 
    ### block 4 ###    
    
    ### manual multi-threading ###


    for j in ce11.tss;do
        for i in `seq 1 200`;do
        mA=`less $output_dir/$j.bin_$i.intersect.bed | wc -l`
        A=`grep -i -o a $output_dir/$j.bin_$i.fasta | wc -l`
        FreqSum=`less $output_dir/$j.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`

        OccuSum_M1=`awk '{sum += $4*$8} END {print sum}' $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M1.bed`

        OccuSum_M2=`awk '{sum += $4*$8} END {print sum}' $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M2.bed`

        echo -e "$j.bin_$i\t${mA}\t${A}\t${FreqSum}\t${OccuSum_M1}\t${OccuSum_M2}" >> $output_dir/tss.mAlevel.occu.txt
        done
    done 

###  3.2 6mA level and nucleosome occupancy ###

### 3. binned TSS enrichment ###


# nohup bash run.tss.optimize.sh >log/tss/run.log 2>&1 &

