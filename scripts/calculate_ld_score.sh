#! /bin/bash
# Calculate LD Scores for Annotations 

input=$1

    echo ${input};
    chr=`echo $input | cut -d . -f 3 <<< $input`
    prefix=`echo $input | cut -d . -f 1 <<< $input`
    echo $input
    python /hpc/packages/minerva-common/ldsc/1.0.0/ldsc/ldsc.py \
    --l2 \
    --bfile /sc/arion/projects/ad-omics/ashvin/ldsc_annotations/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
    --ld-wind-cm 1 \
    --annot ${input}.annot.gz \
    --out ${input}

    zcat ${input}.annot.gz | cut -f 6 > ${input}.annot;
    rm ${input}.annot.gz;
    gzip ${input}.annot;
    
    zcat ${input}.l2.ldscore.gz | cut -f 1,2,3,5 > ${input}.l2.ldscore;
    rm ${input}.l2.ldscore.gz; 
    gzip ${input}.l2.ldscore;
    
    cat ${input}.l2.M_5_50 | cut -f 2 > ${input}.l2.M_5_50.revised;
    rm ${input}.l2.M_5_50;
    mv ${input}.l2.M_5_50.revised ${input}.l2.M_5_50;

    cat ${input}.l2.M | cut -f 2 > ${input}.l2.M.revised;
    rm ${input}.l2.M;
    mv ${input}.l2.M.revised ${input}.l2.M;
    
    rm ${input}.log;
   

