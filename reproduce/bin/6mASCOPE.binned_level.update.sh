#!/bin/sh 

set -eou pipefail


SAMPLE_NAME=Algae
SAMPLE_GENOME_NAME=C.reinhardtii
MAPPED_DIR=/home/user/data2/lit/DATA/project/6mA/reproduce/output/6mASCOPE/Algae/Algae.ccs.fasta.Creinhardtii_281_v5.0.fa.mapped
OUTPUT_DIR=/home/user/data2/lit/DATA/project/6mA/reproduce/output/6mASCOPE/Algae/center
BINNED_BED_DIR=/home/user/data2/lit/project/6mA/reproduce/output/feature/binned/center/C.reinhardtii
# chr start end bin_number
IPD_OUT=/home/user/data2/lit/project/6mA/reproduce/6mASCOPE/Algae.Ipd.out
TMP_DIR=/home/user/data2/lit/DATA/project/6mA/reproduce/output/6mASCOPE/Algae/tmp/center
PREDICT_SCRIPT=/home/user/data2/lit/project/6mA/reproduce/6mASCOPE/code/predict_6mA_level.sh

[ -d $TMP_DIR ] || mkdir $TMP_DIR
[ -d $OUTPUT_DIR ] || mkdir $OUTPUT_DIR
<<'!'
# unmapped reads are already filtered

# some reads FLAG==16 || FLAG==0 but MAPQ==0 ( representative alignment )
awk '{if ($5>0) print $0}' $MAPPED_DIR/$SAMPLE_NAME.ccs.fasta.goi.ref.minimap2.mapped > \
$MAPPED_DIR/$SAMPLE_NAME.ccs.fasta.goi.ref.minimap2.mapped.filter.sam

# filter supplementary alignment   ( FLAG ==2048 )
awk '{if($2==16 || $2==0) print $0}' $MAPPED_DIR/$SAMPLE_NAME.ccs.fasta.goi.ref.minimap2.mapped.filter.sam > \
$MAPPED_DIR/$SAMPLE_NAME.ccs.fasta.goi.ref.minimap2.mapped.filter.2048.sam

# ref and ccs relative position
less $MAPPED_DIR/$SAMPLE_NAME.ccs.fasta.goi.ref.minimap2.mapped.filter.2048.sam | \
awk -v OFS='\t' '{if ($2==16) print $3,$4,$4+length($10)-1,$1,length($10),1,"-";else if($2==0) print $3,$4,$4+length($10)-1,$1,1,length($10),"+"}' | \
sort -k1,1 -k2,2n >$OUTPUT_DIR/$SAMPLE_NAME.paf_like.bed
!

# cp /home/user/data2/lit/DATA/project/6mA/reproduce/output/6mASCOPE/Algae/Algae.paf_like.bed /home/user/data2/lit/DATA/project/6mA/reproduce/output/6mASCOPE/Algae/center/
# ipd.out.A.bed
# less $IPD_OUT | awk -v OFS='\t' '{if($4=="A" && $3==0) print $1,$2-1,$2,$4,$5,"+",$0 ; if($4=="A" && $3==1) print $1,$2-1,$2,$4,$5,"-",$0}' | sort -k1,1 -k2,2n > $OUTPUT_DIR/$SAMPLE_NAME.Ipd.out.A.bed 

<<'!'
time ( 
for i in `seq 1 10`;do
    for j in nucleosome linker;do
        echo "*** intersect first time  ${j}_${i}"
        bedtools intersect -wb -a $BINNED_BED_DIR/$SAMPLE_GENOME_NAME.$j.bin_$i.bed -b $OUTPUT_DIR/$SAMPLE_NAME.paf_like.bed -sorted > $TMP_DIR/$SAMPLE_GENOME_NAME.$j.bin_$i.intersect.bed &
    done
done
wait )
!

for bin in `ls $BINNED_BED_DIR/*bin_*bed`;do
prefix_bed=`basename $bin` && prefix=${prefix_bed%.*}
echo $prefix
done > $OUTPUT_DIR/binned_bed.lst

time (
while read bin;do
    echo "*** intersect first time $bin"
    bedtools intersect -wb -a $BINNED_BED_DIR/$bin.bed -b $OUTPUT_DIR/$SAMPLE_NAME.paf_like.bed -sorted > $TMP_DIR/$bin.intersect.bed &
done < $OUTPUT_DIR/binned_bed.lst
wait )

<<'!'
for i in `seq 1 10`;do
    for j in nucleosome linker;do
        awk -v OFS='\t' '{if ($11=="+") print $8,$9+($2-$6)-1,$10-($7-$3),".",".",$11;if($11=="-") print $8,$10+($7-$3)-1,$9-($2-$6),".",".",$11}' $TMP_DIR/$SAMPLE_GENOME_NAME.$j.bin_$i.intersect.bed | \
        sort -k1,1 -k2,2n > $TMP_DIR/$SAMPLE_GENOME_NAME.$j.bin_$i.intersect.ccs.bed
    done
done
!

less $OUTPUT_DIR/binned_bed.lst | \
while read bin;do
    awk -v OFS='\t' '{if ($11=="+") print $8,$9+($2-$6)-1,$10-($7-$3),".",".",$11;if($11=="-") print $8,$10+($7-$3)-1,$9-($2-$6),".",".",$11}' $TMP_DIR/$bin.intersect.bed | \
    sort -k1,1 -k2,2n > $TMP_DIR/$bin.intersect.ccs.bed
done

<<'!'
time ( 
for i in `seq 1 10`;do
    for j in nucleosome linker;do
        echo "*** intersect second time ${j}_${i}"
        # bedtools intersect -wa -wb -a Algae.Ipd.out.A.bed -b  $SAMPLE_GENOME_NAME.$j.bin_$i.intersect.bed -sorted -s > $SAMPLE_GENOME_NAME.$j.bin_$i.ipd.out.A
        bedtools intersect -wa -wb -a $OUTPUT_DIR/$SAMPLE_NAME.Ipd.out.A.bed -b  $TMP_DIR/$SAMPLE_GENOME_NAME.$j.bin_$i.intersect.ccs.bed -s > $TMP_DIR/$SAMPLE_GENOME_NAME.$j.bin_$i.ipd.out.A.bed &
    done
done
wait )
!

time (
while read bin;do
    echo "*** intersect second time $bin"
    bedtools intersect -wa -wb -a $OUTPUT_DIR/$SAMPLE_NAME.Ipd.out.A.bed -b  $TMP_DIR/$bin.intersect.ccs.bed -s > $TMP_DIR/$bin.ipd.out.A.bed &
done < $OUTPUT_DIR/binned_bed.lst
wait )

<<'!'
for i in `seq 1 10`;do
    for j in nucleosome linker;do
        cut -f 7-16  $TMP_DIR/$SAMPLE_GENOME_NAME.$j.bin_$i.ipd.out.A.bed > $TMP_DIR/$SAMPLE_GENOME_NAME.$j.bin_$i.ipd.out.A
    done
done
!


less $OUTPUT_DIR/binned_bed.lst | \
while read bin;do
    cut -f 7-16  $TMP_DIR/$bin.ipd.out.A.bed > $TMP_DIR/$bin.ipd.out.A
done


cd $TMP_DIR
[ -d merge ] || mkdir -p merge 
for j in linker;do
    cat $SAMPLE_GENOME_NAME.$j.bin_5.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_6.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_5.ipd.out.A
    cat $SAMPLE_GENOME_NAME.$j.bin_4.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_7.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_4.ipd.out.A
    cat $SAMPLE_GENOME_NAME.$j.bin_3.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_8.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_3.ipd.out.A
    cat $SAMPLE_GENOME_NAME.$j.bin_4.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_9.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_2.ipd.out.A
    cat $SAMPLE_GENOME_NAME.$j.bin_1.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_10.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_1.ipd.out.A
done

for j in nucleosome;do
    cat $SAMPLE_GENOME_NAME.$j.bin_5.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_6.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_1.ipd.out.A
    cat $SAMPLE_GENOME_NAME.$j.bin_4.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_7.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_2.ipd.out.A
    cat $SAMPLE_GENOME_NAME.$j.bin_3.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_8.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_3.ipd.out.A
    cat $SAMPLE_GENOME_NAME.$j.bin_4.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_9.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_4.ipd.out.A
    cat $SAMPLE_GENOME_NAME.$j.bin_1.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_10.ipd.out.A > merge/$SAMPLE_GENOME_NAME.$j.bin_5.ipd.out.A
done

[ -f $OUTPUT_DIR/binned.level.txt ] && rm -rf $OUTPUT_DIR/binned.level.txt
[ -f $OUTPUT_DIR/binned.level.merge.txt ] && rm -rf $OUTPUT_DIR/binned.level.merge.txt 

. "/home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh"
conda activate predict_6mASCOPE

<<'!'
# & lead to errors? data disk error-prone
time (
for j in nucleosome linker;do 
    for i in `seq 1 5`;do
        echo "*** bin5 predict ${j}_${i}"
        $PREDICT_SCRIPT merge/$SAMPLE_GENOME_NAME.$j.bin_$i.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_$i >> $OUTPUT_DIR/binned.level.merge.txt 
    done
done
) 
!


time (
less $OUTPUT_DIR/binned_bed.lst | grep -e '[1-5]$' |    \
while read bin;do
    echo "*** bin5 predict $bin"
    $PREDICT_SCRIPT merge/$bin.ipd.out.A $bin >> $OUTPUT_DIR/binned.level.merge.txt 
done
)


<<'!'
time (
for j in nucleosome linker;do
    for i in `seq 1 10`;do
        echo "*** bin10 predict ${j}_${i}"
        $PREDICT_SCRIPT $SAMPLE_GENOME_NAME.$j.bin_$i.ipd.out.A $SAMPLE_GENOME_NAME.$j.bin_$i >> $OUTPUT_DIR/binned.level.txt 
    done
done
)
!


time (
less $OUTPUT_DIR/binned_bed.lst | \
while read bin;do
    echo "*** bin10 predict $bin"
    $PREDICT_SCRIPT $bin.ipd.out.A $bin >> $OUTPUT_DIR/binned.level.txt 
done
)

cd -

echo "done"


# nohup bash /home/user/data2/lit/project/6mA/reproduce/bin/6mASCOPE.binned_level.update.sh > /home/user/data2/lit/project/6mA/reproduce/log/6mASCOPE.binned_level.center.log 2>&1 &