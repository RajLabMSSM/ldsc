#! /bin/bash
#liftover bed files to hg19

#ml util-linux
#input=$1
cd bed_files
ml liftover/09-Jul-2019

#for f in $(ls *.bed); do cut -f 2,3,4,5 $f > decoy_$f; done

for f in $(ls caQTL_*.bed); do liftOver $f /sc/arion/projects/ad-omics/ashvin/hg38ToHg19.over.chain.gz hg19_$f unMapped_$f ; done


#for f in $(ls *.l2.ldscore.gz); do zcat $f | awk -i inplace -F'|' '{print $1,$2,$3,$5}' $f; done
#gzip *ldscore

