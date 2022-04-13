# #!/bin/sh

set -eou pipefail

subreadBam=$1
basename=`basename $subreadBam`
subreadBamName=${basename%%.*}
# output_dir=$2/${subreadBamName}
# log_dir=$3

# [ -d ${output_dir} ] || mkdir -p ${output_dir}
# [ -d ${log_dir} ] || mkdir -p ${log_dir}

time ( 6mASCOPE ccs -i $subreadBam -ob $subreadBamName.ccs.bam -of $subreadBamName.ccs.fasta -p 50 ) > \
ccs.${subreadBamName}.log 2>&1

time ( 6mASCOPE ipd -r $subreadBam -f $subreadBamName.ccs.fasta -o $subreadBamName.Ipd.out -p 50 ) > \
ipd.${subreadBamName}.log 2>&1
