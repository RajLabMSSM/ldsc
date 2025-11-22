#! /bin/bash
#liftover bed files to hg19

#ml util-linux
ml liftover/09-Jul-2019
input_bed=$1
output_bed=$2 
genome_build=$3
#hg19="hg19"
echo $genome_build

if [[ "${genome_build}" == "hg19" ]]
then
   cp $input_bed $output_bed
elif [[ "$genome_build" == "hg38" ]]
then
   liftOver -bedPlus=4 $input_bed data/hg38ToHg19.over.chain.gz $output_bed ${output_bed}_unMapped
   rm ${output_bed}_unMapped
   echo hello!
else 
   echo "Please specify correct hg build - LDSC cannot support non-hg19 or hg38 files."
fi


