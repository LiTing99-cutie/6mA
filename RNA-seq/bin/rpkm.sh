################################################
#File Name: bin/rpkm.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 03 Mar 2022 04:35:26 PM CST
################################################

#!/bin/sh 


### library size (# mapped reads)
[ -d output/rpkm/ ] || mkdir output/rpkm/
output_dir=output/rpkm/
[ -f $output_dir/mappedReadCounts.txt  ] && rm -f $output_dir/mappedReadCounts.txt
mapping_qc_dir=output/mapping_qc

for i in sample;do for j in `seq 1 12`;do echo ${i}${j} ;done ;done |
while read sample;do
    countMap=$( cat $mapping_qc_dir/samtoolsFlagstat.$sample.clean.txt | grep 'mapped (' | sed -E 's@\s\+\s.*@@g' )
    echo -e "$sample\t${countMap}" >> $output_dir/mappedReadCounts.txt
done




#### gene length from gtftools.py
[ -d output/rpkm/gene_length/ ] || mkdir -p output/rpkm/gene_length/
output_dir=output/rpkm/gene_length/

# script dir
gtftools=/home/user/data2/uplee/tools/GTFtools_0.8.0/gtftools.py 

# gunzip -c /home/user/data/lit/database/public/annotation/gtf/ensembl/ensembl.105/Caenorhabditis_elegans.WBcel235.105.gtf.gz > \
# /home/user/data/lit/database/public/annotation/gtf/ensembl/ensembl.105/Caenorhabditis_elegans.WBcel235.105.gtf

gtf=/home/user/data/lit/database/public/annotation/gtf/ensembl/ensembl.105/Caenorhabditis_elegans.WBcel235.105.gtf

gtftools.py -l -c chrI,chrII,chrIII,chrIV,chrV,chrX,chrM $output_dir/ce11.ensembl.105.geneLen.txt $gtf
