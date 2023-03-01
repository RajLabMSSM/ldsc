# LDSC Pipeline

This is the LD Score Regression (LDSC) Pipeline. LDSC quantifies the separate contributions of polygenic effects and various confounding factors, such as population stratification, based on GWAS summary statistics. 

For our purposes, we are interested in producing SNP-based heritability estimates and using stratified LDSC (S-LDSC), and partition this heritability into separate categories. 

This pipeline will create a binary annotation + LD score files given a bed file as an input. It will then assess partitioned heritability using LD Score Regression using the annotations/LD score files created.  

To Run: 

1. First edit the config.yaml file as necessary. The configurations are as follows: 
    - 'GWAS': GWAS summary statistics name (the pipeline uses the QTL-GWAS excel spreadsheet name column for this, ex. Farrell_PSP)
    - 'inFolder': Folder containing all of the bed files you are interested in using 
    - 'prefix': Prefix of your bed file (ex. for file "Microglia_promoters.bed", your prefix would be "Microglia_promoters")
    - 'gwas_inputfile': name for your file that contains paths to all of the LDSC annotations the pipeline will create. Needs to have ".ldcts" as the suffix (ex. if you are using Nott et al. annotations, then you can name your file "Nott.ldcts"
    - 'outfolder': outfolder for all of your results 

2. Activate your snakemake conda environment. 

```
conda activate snakemake
```

3. Do a dry run to make sure you are obtaining the files that you are interested in (both your annotations + LD score files as well as the results for LDSC)

```
snakemake -s Snakefile --configfile config.yaml -npr
```

4. Start your run! Use snakejob_HPC convention. 

```
snakejob_HPC -s Snakefile 
```
