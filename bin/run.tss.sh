

### 3. binned TSS enrichment ###


set -e

output_dir=/home/user/data/lit/project/6mA/iNPS_output
ce11_dir=/home/user/data/lit/database/public/genome/ce11/
ce11_chrSize=$ce11_dir/ce11.fasta.size

### 3.1 bin ###
<<'!'
### bin tss 2kb upstream and downstream ###
    gene_dir=/home/user/data/lit/database/public/annotation/gene_and_gene_predictions/
    grep -v '^#' $gene_dir/ce11.refGene.bed | awk -v OFS='\t' '{print $3,$5-2000,$5+2000}' |  grep -v '-' > $gene_dir/ce11.tss.bed
    bedtools makewindows -b $gene_dir/ce11.tss.bed  -w 10 -i winnum > \
    $output_dir/ce11.tss.binned.bed
### bin tss 2kb upstream and downstream ###
!
### 3.1 bin ###

<<'!'
### like_wiggle to wig ###
    for i in 1 2;do
        for j in chrI chrII chrIII chrIV chrV chrX chrM; do
            awk -v OFS="\t" '{sub(/chromosome/,"chrom"); print $0}'  $output_dir/PC10_LY_MCC-M${i}_$j.like_wig | sed -n '2p' > \
            $output_dir/PC10_LY_MCC-M${i}_$j.like_wig.title

            awk -v OFS="\t" '{print $1,$2}' $output_dir/PC10_LY_MCC-M${i}_$j.like_wig | sed '1,3d' > \
            $output_dir/PC10_LY_MCC-M${i}_$j.like_wig.titleRm.ori

            cat $output_dir/PC10_LY_MCC-M${i}_$j.like_wig.title $output_dir/PC10_LY_MCC-M${i}_$j.like_wig.titleRm.ori > \
            $output_dir/PC10_LY_MCC-M${i}_$j.OriginalProfile.wig

            awk -v OFS="\t" '{print $1,$3}' $output_dir/PC10_LY_MCC-M${i}_$j.like_wig | sed '1,3d' > \
            $output_dir/PC10_LY_MCC-M${i}_$j.like_wig.titleRm.smoo

            cat $output_dir/PC10_LY_MCC-M${i}_$j.like_wig.title $output_dir/PC10_LY_MCC-M${i}_$j.like_wig.titleRm.smoo > \
            $output_dir/PC10_LY_MCC-M${i}_$j.SmoothProfile.wig
        done
    done
### like_wiggle to wig ###
!

### wig to bed; 10-bin bed to single-base bed ###

    # for i in 1 2;do
    #     for profile in OriginalProfile SmoothProfile;do
    #         [ -f $output_dir/PC10_LY_MCC-M${i}.$profile.singlebase.bed ] && rm -rf $output_dir/PC10_LY_MCC-M${i}.$profile.singlebase.bed
    #     done
    # done

    # time ( 
    # for i in 1 2;do
    #     for j in chrI chrII chrIII chrIV chrV chrX chrM; do
    #         for profile in OriginalProfile SmoothProfile;do
    #             echo "wig2bed PC10_LY_MCC-M${i}_$j.$profile"
    #             wig2bed < $output_dir/PC10_LY_MCC-M${i}_$j.$profile.wig | awk -v OFS="\t" '{print $1,$2,$3,$5}' > \
    #             $output_dir/PC10_LY_MCC-M${i}_$j.$profile.bed
    #         done
    #     done
    # done ) &

    # time ( 
    # for i in 1 2;do
    # # for i in 1;do
    #     for profile in OriginalProfile SmoothProfile;do
    #     # for profile in OriginalProfile;do
    #         for j in chrI chrII chrIII chrIV chrV chrX chrM; do
    #         # for j in chrI; do
    #             echo "makewindows PC10_LY_MCC-M${i}_$j.$profile"
    #             bedtools makewindows -b $output_dir/PC10_LY_MCC-M${i}_$j.$profile.bed -w 1 -i src > $output_dir/PC10_LY_MCC-M${i}_$j.$profile.binned.bed 
    #             # head -n 40000 $output_dir/PC10_LY_MCC-M${i}_$j.$profile.bed | bedtools makewindows -b - -w 1 -i src > $output_dir/PC10_LY_MCC-M${i}_$j.$profile.binned.bed
    #         done
    #     done
    # done ) 

    # the server restarted
    # time ( 
    # for i in 2;do
    # # for i in 1;do
    #     for profile in SmoothProfile;do
    #     # for profile in OriginalProfile;do
    #         for j in chrV chrX chrM; do
    #         # for j in chrI; do
    #             echo "makewindows PC10_LY_MCC-M${i}_$j.$profile"
    #             bedtools makewindows -b $output_dir/PC10_LY_MCC-M${i}_$j.$profile.bed -w 1 -i src > $output_dir/PC10_LY_MCC-M${i}_$j.$profile.binned.bed 
    #             # head -n 40000 $output_dir/PC10_LY_MCC-M${i}_$j.$profile.bed | bedtools makewindows -b - -w 1 -i src > $output_dir/PC10_LY_MCC-M${i}_$j.$profile.binned.bed
    #         done
    #     done
    # done ) 

    # for i in 1 2;do
    #     for profile in OriginalProfile SmoothProfile;do
    #         for j in chrI chrII chrIII chrIV chrV chrX chrM; do
    #             cat $output_dir/PC10_LY_MCC-M${i}_$j.$profile.binned.bed >> $output_dir/PC10_LY_MCC-M${i}.$profile.singlebase.bed
    #         done
    #     done
    # done 
### wig to bed; 10-bin bed to single-base bed ###


###  3.2 6mA level and nucleosome occupancy ###
    [ -f tss.mAlevel.occu.txt ] && rm -rf tss.mAlevel.occu.txt
    time ( 
        # for i in 1;do
        for i in `seq 1 400`;do
        for j in ce11.tss;do
            echo "*** process $j.bin_$i"
            less $output_dir/$j.binned.bed | \
            # head -n 40000 | \
            grep $i$ | awk -v OFS='\t' '{print $1,$2,$3}' > $output_dir/$j.bin_$i.bed
            bedtools intersect -wb -a $output_dir/$j.bin_$i.bed -b sig_gt_2_fc0_cov10_for_lt_wt1.bed | \
            awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/$j.bin_$i.intersect.bed 
            bedtools intersect -wb -a $output_dir/$j.bin_$i.bed -b $output_dir/PC10_LY_MCC-M1.SmoothProfile.singlebase.bed | \
            awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M1.bed 
            bedtools intersect -wb -a $output_dir/$j.bin_$i.bed -b $output_dir/PC10_LY_MCC-M2.SmoothProfile.singlebase.bed | \
            awk -v OFS="\t" '{print $4,$5,$6,$7}' > $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M2.bed 
            bedtools getfasta -fi $ce11_dir/ce11.fa -bed $output_dir/$j.bin_$i.bed > $output_dir/$j.bin_$i.fasta
            mA=`less $output_dir/$j.bin_$i.intersect.bed | wc -l`
            A=`grep -i -o a $output_dir/$j.bin_$i.fasta | wc -l`
            FreqSum=`less $output_dir/$j.bin_$i.intersect.bed | awk '{sum += $4} END {print sum}'`
            OccuSum_M1=`less $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M1.bed | awk '{sum += $4} END {print sum}'`
            OccuSum_M2=`less $output_dir/$j.bin_$i.intersect.occu.PC10_LY_MCC-M2.bed | awk '{sum += $4} END {print sum}'`
            echo -e "$j.bin_$i\t${mA}\t${A}\t${FreqSum}\t${OccuSum_M1}\t${OccuSum_M2}" >> tss.mAlevel.occu.txt
        done
    done ) 

###  3.2 6mA level and nucleosome occupancy ###

### 3. binned TSS enrichment ###


# nohup bash run.tss.sh >log/feature.2.log 2>&1 &
# nohup bash run.tss.sh >log/feature.2.1.log 2>&1 &
# nohup bash run.tss.sh >log/feature.2.2.log 2>&1 &
