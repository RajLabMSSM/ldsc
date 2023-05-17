#! /bin/bash
#liftover bed files to hg19

#ml util-linux
# ml liftover/09-Jul-2019
input_bed=$1
output_bed=$2 
genome_build=$3
hg19="hg19"
echo $genome_build

if [[ "$genome_build"=="${hg19}" ]];
then
   cp $input_bed $output_bed
else
   liftOver $input_bed data/hg38ToHg19.over.chain.gz $output_bed ${output_bed}_unMapped
   rm ${output_bed}_unMapped
fi 
