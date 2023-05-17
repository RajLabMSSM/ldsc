import sys 
import os.path
import glob
import pandas as pd 

gwas = config["GWAS"]
chromosome_numbers = list(range(1, 23))
infolder = config["inFolder"]
annot_prefix = config["prefix"]
outFolder = config["outFolder"]
# revised_bed_file = ["revised_" + b for b in bed_file]

ldcts_prefix = config["gwas_inputfile"]
gwas_df = pd.read_excel("/sc/arion/projects/ad-omics/data/references//GWAS/GWAS-QTL_data_dictionary.xlsx", sheet_name = 2).set_index("dataset")
genomeBuild = config['annotation_build']
category = ''.join(ldcts_prefix).split('.')[0]

# ruleorder: create_annotation > split_by_chr > calculate_ld_score > modify_cts_annotations > new_ldcts_file > format_sumstats > munge_sumstats

rule all: 
   input:
      # expand("{bedfile_input}", bedfile_input = bed_file_input)
      # expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.ldscore.gz", annotprefix = annot_prefix, chr_nums = chromosome_numbers, outfolder = outFolder), 
      # expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M", annotprefix = annot_prefix, chr_nums = chromosome_numbers, outfolder = outFolder),
      # expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M_5_50", annotprefix = annot_prefix, chr_nums = chromosome_numbers, outfolder = outFolder)
      # expand("{ldctsprefix}", ldctsprefix = ldcts_prefix)
      # expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.annot.gz", outfolder = outFolder, annotprefix = annot_prefix, chr_nums = chromosome_numbers)
      expand("{outfolder}/enrichment_results/{GWAS}.cell_type_results.txt", outfolder=outFolder, GWAS=gwas)
      # expand("bed_files/annotation_files/{bedfile}.a.{chr_nums}.annot.gz.revised.annot.gz", chr_nums = chromosome_numbers, bedfile = bed_file),
      # expand("bed_files/annotation_files/{bedfile}.a.{chr_nums}.l2.ldscore.gz.revised.l2.ldscore.gz", chr_nums = chromosome_numbers, bedfile = bed_file),
      # expand("bed_files/annotation_files/{bedfile}.a.{chr_nums}.l2.M.revised.l2.M", chr_nums = chromosome_numbers, bedfile = bed_file),
      # expand("bed_files/annotation_files/{bedfile}.a.{chr_nums}.l2.M_5_50.revised.l2.M_5_50", chr_nums = chromosome_numbers, bedfile = bed_file) 

rule liftOver_bed:
   input:
      expand("{inFolder}/{{annotprefix}}.bed", inFolder = infolder, annotprefix = annot_prefix)
   output:
      expand("annotation_database/{c}/input_bed_files/{{annotprefix}}.bed", c = category, annotprefix = annot_prefix)
   shell:
      "ml liftover/09-Jul-2019;"
      "sh scripts/liftover_script.sh {infolder}/{wildcards.annotprefix}.bed annotation_database/{category}/input_bed_files/{wildcards.annotprefix}.bed {genomeBuild}"
      
      
rule create_annotation:
   input: 
      expand("annotation_database/{c}/input_bed_files/{{annotprefix}}.bed", c = category, annotprefix = annot_prefix)
   output: 
      expand("annotation_database/{c}/annotation_files/{{annotprefix}}.annot", c = category, annotprefix = annot_prefix)
   params: 
      script1="scripts/create_annotation.R"
   shell: 
      "ml R/4.0.3;"
      "Rscript {params.script1} -b {input} -o {output};"

rule split_by_chr: 
   input: 
      expand("annotation_database/{c}/annotation_files/{{annotprefix}}.annot", c = category, annotprefix = annot_prefix)
   output: 
      expand("annotation_database/{c}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.annot.gz", c = category, chr_nums = chromosome_numbers, annotprefix = annot_prefix)
   params:
      c = category
   shell:
      "python scripts/split_by_chr.py --annotation {input} --output annotation_database/{params.c}/annotation_files/{wildcards.annotprefix} --chrom {wildcards.chr_nums}"


rule calculate_ld_score: 
   input:
      expand("annotation_database/{c}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.annot.gz", c = category, chr_nums = chromosome_numbers, annotprefix = annot_prefix)
   output:
      expand("annotation_database/{c}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.l2.ldscore.gz", c = category, chr_nums = chromosome_numbers, outfolder = outFolder), 
      expand("annotation_database/{c}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.l2.M", c = category, chr_nums = chromosome_numbers, annotprefix = annot_prefix),
      expand("annotation_database/{c}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.l2.M_5_50", c = category, chr_nums = chromosome_numbers, outfolder = outFolder, annotprefix = annot_prefix)
   params:
      chr_nums = chromosome_numbers,
      c = category
   shell:
      "ml ldsc/1.0.0;"
      "sh scripts/calculate_ld_score.sh annotation_database/{params.c}/annotation_files/{wildcards.annotprefix}.a.{wildcards.chr_nums};"

rule new_ldcts_file:   
   input:
      expand("annotation_database/{c}/annotation_files/{annotprefix}.a.{chr_nums}.annot.gz", c = category, chr_nums = chromosome_numbers, annotprefix = annot_prefix),
      expand("annotation_database/{c}/annotation_files/{annotprefix}.a.{chr_nums}.l2.ldscore.gz", c = category, chr_nums = chromosome_numbers, outfolder = outFolder, annotprefix = annot_prefix), 
      expand("annotation_database/{c}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M", c = category, chr_nums = chromosome_numbers, annotprefix = annot_prefix),
      expand("annotation_database/{c}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M_5_50", c = category, chr_nums = chromosome_numbers, outfolder = outFolder, annotprefix = annot_prefix)   
   output:
      expand("{ldctsprefix}", ldctsprefix = ldcts_prefix, outfolder = outFolder)
   params: 
      ldctsprefix = ldcts_prefix,
      outfolder = outFolder,
      script1 = "scripts/create_ldcts.R",
      annotprefix = annot_prefix,
      joined_bar=lambda w, input: ",".join(annot_prefix),  # ', input' was added
      c = category
   shell:
      # "echo ${params.annotprefix}"
      "ml R/4.0.3;"
      "Rscript {params.script1} -i annotation_database/{params.c}/annotation_files/ -b {params.joined_bar} -o {params.ldctsprefix};"

rule format_sumstats: 
   output:
      expand("{outfolder}/formatted_ldsc_gwas/{GWAS}.sumstats.gz", outfolder = outFolder, GWAS=gwas)
   params:
      out_folder = outFolder
   run:
      import numpy as np
      gwas_dict = pd.read_excel("/sc/arion/projects/ad-omics/data/references//GWAS/GWAS-QTL_data_dictionary.xlsx", sheet_name = 2)
      def format_gwas(gwas_d): 
         new_gwas_dict = gwas_dict.loc[gwas_dict['dataset'] == gwas_d]
         print(new_gwas_dict['full_processed_path'].values[0])
         gwas_df = pd.read_csv(new_gwas_dict['full_processed_path'].values[0], sep = '\t')
         column_names = [new_gwas_dict['full_snp'].values[0], 
                    new_gwas_dict['full_A1'].values[0], 
                    new_gwas_dict['full_A2'].values[0], 
                    new_gwas_dict['full_p'].values[0],
                    new_gwas_dict['full_effect'].values[0],
                    new_gwas_dict['full_se'].values[0]]
         
         gwas_shortened = gwas_df[column_names]
         print(new_gwas_dict['full_effect'].str.contains("OR"))
	 print(gwas_shortened.iloc[:,4], gwas_shortened.iloc[:,5])
         if (new_gwas_dict['full_effect'].str.contains("OR").values[0] == True):
            gwas_shortened['Z'] = gwas_shortened.iloc[:,4] - 1  
            # gwas_shortened['Z'] = np.log10(gwas_shortened.iloc[:,4]) / gwas_shortened.iloc[:,5]
         else: 
            gwas_shortened['Z'] = gwas_shortened.iloc[:,4] / gwas_shortened.iloc[:,5]
  
         gwas_shortened = gwas_shortened.iloc[:, [0, 1, 2, 3, 6]]
         gwas_shortened.columns = ['ID', 'Allele1', 'Allele2', 'P.value', 'Z_score']
         gwas_shortened.dropna(inplace=True)
         print(gwas_shortened.head())
         return(gwas_shortened)
      print("hello!")
      print(gwas)
      for g in gwas: 
         formatted_gwas = format_gwas(g)
	 print("goodbye!")
         formatted_gwas.to_csv(''.join(params.out_folder) + "/formatted_ldsc_gwas/" + g + '.sumstats.gz', sep = '\t', index = None)

rule munge_sumstats:
   input: 
      expand("{{outfolder}}/formatted_ldsc_gwas/{{GWAS}}.sumstats.gz", outfolder=outFolder, GWAS=gwas)
   output:
      expand("{{outfolder}}/munged_ldsc_gwas/{{GWAS}}_munged.sumstats.gz", outfolder=outFolder, GWAS=gwas)
   params:
      sample_size = lambda wildcards: gwas_df.loc[wildcards.GWAS]['N'],
   shell:
      "ml ldsc/1.0.0;"
      "python /hpc/packages/minerva-common/ldsc/1.0.0/ldsc/munge_sumstats.py --sumstats {wildcards.outfolder}/formatted_ldsc_gwas/{wildcards.GWAS}.sumstats.gz --merge-alleles  data/w_hm3.snplist --out {wildcards.outfolder}/munged_ldsc_gwas/{wildcards.GWAS}_munged --N {params.sample_size} --snp ID --a1 Allele1 --a2 Allele2 --p P.value --signed-sumstats Z_score,0"

rule run_ldsc:
   input:
      expand("{ldctsprefix}", ldctsprefix = ldcts_prefix),
      expand("{outfolder}/munged_ldsc_gwas/{{GWAS}}_munged.sumstats.gz", outfolder=outFolder, GWAS=gwas)
   output:
      expand("{outfolder}/enrichment_results/{{GWAS}}.cell_type_results.txt", outfolder=outFolder, GWAS=gwas)
   shell:
      "ml ldsc/1.0.0;"
      "python /hpc/packages/minerva-common/ldsc/1.0.0/ldsc/ldsc.py --h2-cts {outFolder}/munged_ldsc_gwas/{wildcards.GWAS}_munged.sumstats.gz --ref-ld-chr /sc/arion/projects/ad-omics/ashvin/ldsc_annotations/1000G_EUR_Phase3_baseline/baseline. --out {outFolder}/enrichment_results/{wildcards.GWAS} --ref-ld-chr-cts {ldcts_prefix} --w-ld-chr /sc/arion/projects/ad-omics/ashvin/ldsc_annotations/weights_hm3_no_hla/weights."
