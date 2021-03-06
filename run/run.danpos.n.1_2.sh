################################################
#File Name: /home/user/data2/lit/project/6mA/run.danpos.n.1_2.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 16 Mar 2022 03:53:28 PM CST
################################################

#!/bin/sh 


set -e

<<'!'
###  feature analysis ###

    [ -d /home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2 ] || mkdir -p /home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2
    [ -d log/danpos/normed/method_1_2 ] || mkdir -p log/danpos/normed/method_1_2
    output_dir=/home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2
    log_dir=log/danpos/normed/method_1_2
    ce11_dir=/home/user/data/lit/database/public/genome/ce11/pengq/
    ce11_chrSize=$ce11_dir/Caenorhabditis_elegans.WBcel235.dna.fa.size

    ### nucleosome bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.nucleosome.bed"
        sed '1d' danpos/normed/${sample}/pooled/out.uniq.rmdup.nameSorted.smooth.positions.xls | \
        # nucleosomes are all kept
        # sorted for bedtools complement 
        awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $output_dir/${sample}.nucleosome.bed

        nucleosome_number=`less $output_dir/${sample}.nucleosome.bed | wc -l`
        echo "${sample} nucleosome_number: ${nucleosome_number}"
    done )
    ### nucleosome bed ###
    
    ### linker bed ###
    time (
    for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
        echo "*** generating ${sample}.linker.bed"
        # do not remove the first linker and the last linker
        bedtools complement  -i $output_dir/${sample}.nucleosome.bed -g $ce11_chrSize > \
        $output_dir/${sample}.linker.bed

        linker_number=`less $output_dir/${sample}.linker.bed | wc -l`
        echo "${sample} linker_number: ${linker_number}"
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
        done 
        wait )

    ### 2.1 bin ###

    ### 2.2 6mA level ###
        [ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt
        time (
        for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
            for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
                for j in nucleosome linker;do
                    for i in `seq 1 10`;do
                        list_tmp=${mA_list%.*} && list=${list_tmp##*_} 
                        echo "*** process ${sample}.$j.${list}.bin_$i"
                        less $output_dir/${sample}.$j.binned.bed | \
                                grep $i$ > \
                                    $output_dir/${sample}.$j.bin_$i.bed
                                
                        less $output_dir/${sample}.$j.bin_$i.bed | awk -v OFS='\t' '{if($1=="MtDNA") print "chrM",$2,$3,$4 ; else print "chr"$1,$2,$3,$4}' > \
                                    $output_dir/${sample}.$j.bin_$i.chr.bed

                        bedtools intersect -wa -wb -a $output_dir/${sample}.$j.bin_$i.chr.bed -b ${mA_list} > $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed &

                        bedtools nuc -fi $ce11_dir/Caenorhabditis_elegans.WBcel235.dna.fa -bed $output_dir/${sample}.$j.bin_$i.bed > \
                        $output_dir/${sample}.$j.bin_$i.txt &

                    done
                    wait
                done 
            done
        done ) 


        for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
            for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
                for j in nucleosome linker;do
                    for i in `seq 1 10`;do
                        list_tmp=${mA_list%.*} && list=${list_tmp##*_}
                        mA=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | wc -l`
                        A=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' `
                        T=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' `
                        AT=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' `
                        GC=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' `
                        FreqSum=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $8} END {print sum}'`
                        echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${T}\t${FreqSum}\t${AT}\t${GC}" >> $output_dir/binned.level.txt
                    done
                done
            done
        done
    ### 2.2 6mA level ###

    ### 2. binned linker or nucleosome enrichment ###

### feature analysis ###

# 2022-03-16
# [ -d log/danpos/normed/method_1_2 ] || mkdir -p log/danpos/normed/method_1_2
# nohup bash run.danpos.n.1_2.sh > log/danpos/normed/method_1_2/run.log 2>&1 &
!

<<'!'
awk -v OFS='\t' '{print $1,$2,$3,$4,NR}' /home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2/PC10_LY_MCC-M1.nucleosome.bin_1.chr.bed > \
/home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2/PC10_LY_MCC-M1.nucleosome.bin_1.chr.num.bed
cut -f 1-3 /home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2/PC10_LY_MCC-M1.nucleosome.wt1.bin_1.intersect.bed | \
bedtools intersect -wa -wb -a /home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2/PC10_LY_MCC-M1.nucleosome.bin_1.chr.num.bed \
-b -
!

output_dir=/home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2
[ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt
for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
        for j in nucleosome linker;do
            for i in `seq 1 10`;do
                list_tmp=${mA_list%.*} && list=${list_tmp##*_}
                mA=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | wc -l`
                A=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' `
                C=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' `
                G=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' `
                T=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' `
                base_number=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $13} END {print sum}' `
                AT=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' `
                GC=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' `
                FreqSum=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $8} END {print sum}'`
                echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${C}\t${G}\t${T}\t${base_number}\t${AT}\t${GC}\t${FreqSum}" >> $output_dir/binned.level.txt
            done
        done
        
    done
done

# nohup bash run.danpos.n.1_2.sh >> log/danpos/normed/method_1_2/run.log 2>&1 &


# add in 20220412

output_dir=/home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2
mkdir -p $output_dir/merge   

for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
        for j in nucleosome;do
            list_tmp=${mA_list%.*} && list=${list_tmp##*_}
            cat $output_dir/${sample}.$j.${list}.bin_5.intersect.bed $output_dir/${sample}.$j.${list}.bin_6.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_1.intersect.bed
            cat $output_dir/${sample}.$j.${list}.bin_4.intersect.bed $output_dir/${sample}.$j.${list}.bin_7.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_2.intersect.bed
            cat $output_dir/${sample}.$j.${list}.bin_3.intersect.bed $output_dir/${sample}.$j.${list}.bin_8.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_3.intersect.bed
            cat $output_dir/${sample}.$j.${list}.bin_2.intersect.bed $output_dir/${sample}.$j.${list}.bin_9.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_4.intersect.bed
            cat $output_dir/${sample}.$j.${list}.bin_1.intersect.bed $output_dir/${sample}.$j.${list}.bin_10.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_5.intersect.bed
            cat $output_dir/${sample}.$j.bin_5.txt $output_dir/${sample}.$j.bin_6.txt > $output_dir/merge/${sample}.$j.bin_1.txt
            cat $output_dir/${sample}.$j.bin_4.txt $output_dir/${sample}.$j.bin_7.txt > $output_dir/merge/${sample}.$j.bin_2.txt
            cat $output_dir/${sample}.$j.bin_3.txt $output_dir/${sample}.$j.bin_8.txt > $output_dir/merge/${sample}.$j.bin_3.txt
            cat $output_dir/${sample}.$j.bin_2.txt $output_dir/${sample}.$j.bin_9.txt > $output_dir/merge/${sample}.$j.bin_4.txt
            cat $output_dir/${sample}.$j.bin_1.txt $output_dir/${sample}.$j.bin_10.txt > $output_dir/merge/${sample}.$j.bin_5.txt
        done    
    done
done

for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
        for j in linker;do
            list_tmp=${mA_list%.*} && list=${list_tmp##*_}
            cat $output_dir/${sample}.$j.${list}.bin_5.intersect.bed $output_dir/${sample}.$j.${list}.bin_6.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_5.intersect.bed
            cat $output_dir/${sample}.$j.${list}.bin_4.intersect.bed $output_dir/${sample}.$j.${list}.bin_7.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_4.intersect.bed
            cat $output_dir/${sample}.$j.${list}.bin_3.intersect.bed $output_dir/${sample}.$j.${list}.bin_8.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_3.intersect.bed
            cat $output_dir/${sample}.$j.${list}.bin_2.intersect.bed $output_dir/${sample}.$j.${list}.bin_9.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_2.intersect.bed
            cat $output_dir/${sample}.$j.${list}.bin_1.intersect.bed $output_dir/${sample}.$j.${list}.bin_10.intersect.bed > \
            $output_dir/merge/${sample}.$j.${list}.bin_1.intersect.bed
            cat $output_dir/${sample}.$j.bin_5.txt $output_dir/${sample}.$j.bin_6.txt > $output_dir/merge/${sample}.$j.bin_5.txt
            cat $output_dir/${sample}.$j.bin_4.txt $output_dir/${sample}.$j.bin_7.txt > $output_dir/merge/${sample}.$j.bin_4.txt
            cat $output_dir/${sample}.$j.bin_3.txt $output_dir/${sample}.$j.bin_8.txt > $output_dir/merge/${sample}.$j.bin_3.txt
            cat $output_dir/${sample}.$j.bin_2.txt $output_dir/${sample}.$j.bin_9.txt > $output_dir/merge/${sample}.$j.bin_2.txt
            cat $output_dir/${sample}.$j.bin_1.txt $output_dir/${sample}.$j.bin_10.txt > $output_dir/merge/${sample}.$j.bin_1.txt
        done    
    done
done

output_dir=/home/user/data/lit/project/6mA/feature/danpos/normed/method_1_2/merge
[ -f $output_dir/binned.level.txt ] && rm -rf $output_dir/binned.level.txt
for sample in $(for i in `seq 1 2`;do echo PC10_LY_MCC-M$i;done);do
    for mA_list in sig_gt_2_fc0_cov10_for_lt_wt1.bed sig_gt_2_fc0_cov10_for_lt_wt2.bed sig_gt_2_fc0_cov10_for_lt_wt1-wt2.bed;do
        for j in nucleosome linker;do
            for i in `seq 1 5`;do
                list_tmp=${mA_list%.*} && list=${list_tmp##*_}
                mA=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | wc -l`
                A=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $7} END {print sum}' `
                C=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $8} END {print sum}' `
                G=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $9} END {print sum}' `
                T=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $10} END {print sum}' `
                base_number=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $13} END {print sum}' `
                AT=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $5} END {print sum/NR}' `
                GC=` sed '1d' $output_dir/${sample}.$j.bin_$i.txt | awk '{sum += $6} END {print sum/NR}' `
                FreqSum=`less $output_dir/${sample}.$j.${list}.bin_$i.intersect.bed | awk '{sum += $8} END {print sum}'`
                echo -e "${sample}\t$j\t${list}\tbin_$i\t${mA}\t${A}\t${C}\t${G}\t${T}\t${base_number}\t${AT}\t${GC}\t${FreqSum}" >> $output_dir/binned.level.txt
            done
        done 
    done
done