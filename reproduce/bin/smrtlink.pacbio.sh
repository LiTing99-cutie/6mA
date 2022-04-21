#!/bin/sh

set -eou pipefail
# conda create -n pbalign -c bioconda pbalign

. "/home/user/BGM/lit/anaconda3/etc/profile.d/conda.sh"
conda activate pbalign

subreadBam=$1
basename=`basename $subreadBam`
subreadBamName=${basename%%.*}
ref=$2
ccs_fasta=$3
output_dir=$4/${subreadBamName}
log_dir=$5

[ -d ${output_dir} ] || mkdir -p ${output_dir}
[ -d ${log_dir} ] || mkdir -p ${log_dir}

[ -d ${output_dir}/pbmm2/ccs ] || mkdir -p ${output_dir}/pbmm2/ccs
[ -d ${output_dir}/pbmm2/ref ] || mkdir -p ${output_dir}/pbmm2/ref

[ -d ${log_dir}/pbmm2 ] || mkdir -p ${log_dir}/pbmm2

# pbmm2 1.8.0 (commit v1.7.0-9-g3c16a4d)
# pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml|fa|fq|gz|fofn> [out.aligned.bam|xml]

#   ref.fa|xml|mmi             STR    Reference FASTA, ReferenceSet XML, or
#                                     Reference Index
#   in.bam|xml|fa|fq|gz|fofn   STR    Input BAM, DataSet XML, FASTA, or FASTQ
#   out.aligned.bam|xml        STR    Output BAM or DataSet XML

#  -j,--num-threads           INT    Number of threads to use, 0 means
#                                     autodetection. [0]

# ipdSummary 3.0
# ipdSummary aligned.bam --reference ref.fasta m6A,m4C --gff basemods.gff \
# --csv_h5 kinetics.h5

# time ( pbmm2 align -j 50 \
# $ref \
# $subreadBam \
# ${output_dir}/pbmm2/${subreadBamName}.ref.bam ) \
# > ${log_dir}/pbmm2/pbmm2.${subreadBamName}.ref.log 2>&1 

# nohup pbindex /home/user/data2/lit/project/6mA/reproduce/output/Pacbio/Algae/pbmm2/Algae.ref.bam &

time ( ipdSummary --numWorkers 50 \
${output_dir}/pbmm2/${subreadBamName}.ref.bam \
--reference $ref \
--identify m6A,m4C \
--methylFraction \
--gff ${output_dir}/pbmm2/ref/${subreadBamName}.gff \
--csv ${output_dir}/pbmm2/ref/${subreadBamName}.csv \
--bigwig ${output_dir}/pbmm2/ref/${subreadBamName}.bigwig ) \
> ${log_dir}/pbmm2/ipdSummary.${subreadBamName}.ref.log 2>&1 

# time ( pbmm2 align -j 50 \
# $ccs_fasta \
# $subreadBam \
# ${output_dir}/pbmm2/${subreadBamName}.ccs.bam ) \
# > ${log_dir}/pbmm2/pbmm2.${subreadBamName}.ccs.log 2>&1 

# pbindex /home/user/data2/lit/project/6mA/reproduce/output/Pacbio/Algae/pbmm2/Algae.ccs.bam

# time ( ipdSummary --numWorkers 50 \
# ${output_dir}/pbmm2/${subreadBamName}.ccs.bam \
# --reference $ccs_fasta \
# --identify m6A,m4C \
# --methylFraction \
# --gff ${output_dir}/pbmm2/ccs/${subreadBamName}.gff \
# --csv ${output_dir}/pbmm2/ccs/${subreadBamName}.csv \
# --bigwig ${output_dir}/pbmm2/ccs/${subreadBamName}.bigwig ) \
# > ${log_dir}/pbmm2/ipdSummary.${subreadBamName}.ccs.log 2>&1 