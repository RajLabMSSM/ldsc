import sys 
import os.path
import glob
import pandas as pd 

gwas = config["GWAS"]
chromosome_numbers = list(range(1, 23))
infolder = config["inFolder"]
annot_prefix = config["prefix"]
outFolder = config["outfolder"]
# revised_bed_file = ["revised_" + b for b in bed_file]

ldcts_prefix = config["gwas_inputfile"]
gwas_df = pd.read_excel("/sc/arion/projects/ad-omics/data/references//GWAS/GWAS-QTL_data_dictionary.xlsx", sheet_name = 2)

GWAS_sample_size = [gwas_df.N[gwas_df['dataset'] == g].values[0] for g in gwas]

# ruleorder: create_annotation > split_by_chr > calculate_ld_score > modify_cts_annotations > new_ldcts_file > format_sumstats > munge_sumstats

rule all: 
   input:
      # expand("{bedfile_input}", bedfile_input = bed_file_input)
      # expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.ldscore.gz", annotprefix = annot_prefix, chr_nums = chromosome_numbers, outfolder = outFolder), 
      # expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M", annotprefix = annot_prefix, chr_nums = chromosome_numbers, outfolder = outFolder),
      # expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M_5_50", annotprefix = annot_prefix, chr_nums = chromosome_numbers, outfolder = outFolder)
      # expand("{ldctsprefix}", ldctsprefix = ldcts_prefix)
      # expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.annot.gz", outfolder = outFolder, annotprefix = annot_prefix, chr_nums = chromosome_numbers)
      expand("results/{GWAS}.cell_type_results.txt", GWAS=gwas)
      # expand("bed_files/annotation_files/{bedfile}.a.{chr_nums}.annot.gz.revised.annot.gz", chr_nums = chromosome_numbers, bedfile = bed_file),
      # expand("bed_files/annotation_files/{bedfile}.a.{chr_nums}.l2.ldscore.gz.revised.l2.ldscore.gz", chr_nums = chromosome_numbers, bedfile = bed_file),
      # expand("bed_files/annotation_files/{bedfile}.a.{chr_nums}.l2.M.revised.l2.M", chr_nums = chromosome_numbers, bedfile = bed_file),
      # expand("bed_files/annotation_files/{bedfile}.a.{chr_nums}.l2.M_5_50.revised.l2.M_5_50", chr_nums = chromosome_numbers, bedfile = bed_file) 

rule create_annotation:
   input: 
      expand("{inFolder}/{{annotprefix}}.bed", inFolder = infolder, annotprefix = annot_prefix)
   output: 
      expand("{outfolder}/{{annotprefix}}.annot", outfolder = outFolder, annotprefix = annot_prefix)
   params: 
      script1="create_annotation.R"
   shell: 
      "ml R/4.0.3;"
      "Rscript {params.script1} -b {input} -o {output};"


rule split_by_chr: 
   input: 
      expand("{outfolder}/{{annotprefix}}.annot", outfolder = outFolder, annotprefix = annot_prefix)
   output: 
      expand("{outfolder}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.annot.gz", outfolder = outFolder, chr_nums = chromosome_numbers, annotprefix = annot_prefix)
   params:
      outfolder = outFolder
   shell:
      # "mkdir {params.outfolder}/annotation_files;"
      "python split_by_chr.py --annotation {input} --output {params.outfolder}/annotation_files/{wildcards.annotprefix} --chrom {wildcards.chr_nums}"


rule calculate_ld_score: 
   input:
      expand("{outfolder}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.annot.gz", outfolder = outFolder, chr_nums = chromosome_numbers, annotprefix = annot_prefix)
   output:
      expand("{outfolder}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.l2.ldscore.gz", chr_nums = chromosome_numbers, outfolder = outFolder), 
      expand("{outfolder}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.l2.M", chr_nums = chromosome_numbers, outfolder = outFolder),
      expand("{outfolder}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.l2.M_5_50", chr_nums = chromosome_numbers, outfolder = outFolder),
      # expand("{outfolder}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.log", outfolder = outFolder, annotprefix = annot_prefix, chr_nums = chromosome_numbers)
   params:
      chr_nums = chromosome_numbers,
      outfolder = outFolder
   shell:
      "ml ldsc/1.0.0;"
      "sh calculate_ld_score.sh {params.outfolder}/annotation_files/{wildcards.annotprefix}.a.{wildcards.chr_nums};"
      # "sh modify_ld_score.sh {params.outfolder}/annotation_files/{wildcards.annotprefix}.a.{wildcards.chr_nums}"

   # run:
   #    for f in {params.outfolder}/annotation_files/*.annot.gz; do zcat f | cut -f 1  
   #    annot_files = glob.glob("{params.outfolder}/annotation_files/*.annot.gz")
   #    ldscore_files = glob.glob("{params.outfolder}/annotation_files/*.ldscore.gz")
   #    M_5_50_files = glob.glob("{params.outfolder}/annotation_files/*.l2.M_5_50")
   #    M_files = glob.glob("{params.outfolder}/annotation_files/*.l2.M")
   #    for annot_f in annot_files:
   #       print(annot_f)
   #       annot_frame = pd.read_csv(annot_f, sep = '\t')
   #       revised_annot_frame = annot_frame.iloc[:,5]
   #       revised_annot_frame.to_csv(annot_f + '.revised.annot.gz', compression='gzip', sep='\t', index=None)
   #    for ldscore_f in ldscore_files:
   #       print(ldscore_f)
   #       ldscore_frame = pd.read_csv(ldscore_f, sep = '\t')
   #       revised_ldscore_frame = ldscore_frame.iloc[:,[0,1,2,4]]
   #       revised_ldscore_frame.to_csv(ldscore_f + '.revised.l2.ldscore.gz', compression='gzip', sep='\t', index=None)
   #    for M_5_50 in M_5_50_files:
   #       print(M_5_50)
   #       M_5_50_frame = pd.read_csv(M_5_50, sep = '\t')
   #       revised_M_5_50_frame = M_5_50_frame.iloc[:,1]
   #       revised_M_5_50_frame.to_csv(M_5_50 + '.revised.l2.M_5_50', sep='\t', index=None)
   #    for M in M_files:
   #       print(M)
   #       M_frame = pd.read_csv(M, sep = '\t')
   #       revised_M_frame = M_frame.iloc[:,1]
   #       revised_M_frame.to_csv(M + '.revised.l2.M', sep='\t', index=None)

rule new_ldcts_file:   
   input:
      expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.annot.gz", annotprefix = annot_prefix, outfolder = outFolder, chr_nums = chromosome_numbers),
      expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.ldscore.gz", annotprefix = annot_prefix, outfolder = outFolder, chr_nums = chromosome_numbers),
      expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M", annotprefix = annot_prefix, outfolder = outFolder, chr_nums = chromosome_numbers),
      expand("{outfolder}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M_5_50", annotprefix = annot_prefix, outfolder = outFolder, chr_nums = chromosome_numbers),    
   output:
      expand("{ldctsprefix}", ldctsprefix = ldcts_prefix, annotprefix = annot_prefix)
   params: 
      ldctsprefix = ldcts_prefix,
      outfolder = outFolder,
      script1 = "create_ldcts.R",
      annotprefix = annot_prefix,
      joined_bar=lambda w, input: ",".join(annot_prefix),  # ', input' was added
   shell:
      # "echo ${params.annotprefix}"
      "ml R/4.0.3;"
      "Rscript {params.script1} -i {params.outfolder}/annotation_files/ -b {params.joined_bar} -o {params.ldctsprefix};"
   # run:
   #    def create_ldcts_file(bed_file_list):
   #       print(os.path.basename(bed_file_list))
   #       bed_file_list = bed_file_list.str.split('/').str[0]
   #       print(bed_file_list)
   #       data = pd.DataFrame(bed_file_list)
   #       data['path'] = '[params.annotprefix]/annotation_files/' + bed_file_list + '.a.'
   #       return data
   #    b = annot_prefix
   #    ldcts_file = create_ldcts_file(pd.Series(annot_prefix))
   #    print(''.join(ldcts_prefix))
   #    ldcts_file.to_csv(''.join(ldcts_prefix), sep = '\t', index = None, header = None)


rule format_sumstats: 
   input: 
      expand("{ldctsprefix}", ldctsprefix = ldcts_prefix)
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
         formatted_gwas.to_csv(''.join(params.out_folder) + "formatted_ldsc_gwas/" + g + '.sumstats.gz', sep = '\t', index = None)

rule munge_sumstats:
   input: 
      expand("{{outfolder}}/formatted_ldsc_gwas/{{GWAS}}.sumstats.gz", outfolder=outFolder, GWAS=gwas)
   output:
      expand("{{outfolder}}/munged_ldsc_gwas/{{GWAS}}_munged.sumstats.gz", outfolder=outFolder, GWAS=gwas)
   shell:
      "ml ldsc/1.0.0;"
      "python /hpc/packages/minerva-common/ldsc/1.0.0/ldsc/munge_sumstats.py --sumstats {wildcards.outfolder}/formatted_ldsc_gwas/{wildcards.GWAS}.sumstats.gz --merge-alleles  w_hm3.snplist --out {wildcards.outfolder}/munged_ldsc_gwas/{wildcards.GWAS}_munged --N {GWAS_sample_size} --snp ID --a1 Allele1 --a2 Allele2 --p P.value --signed-sumstats Z_score,0"

rule run_ldsc:
   input:
      expand("{outfolder}/munged_ldsc_gwas/{GWAS}_munged.sumstats.gz", outfolder=outFolder, GWAS=gwas)
   output:
      expand("results/{GWAS}.cell_type_results.txt", GWAS=gwas)
   shell:
      "ml ldsc/1.0.0;"
      "python /hpc/packages/minerva-common/ldsc/1.0.0/ldsc/ldsc.py --h2-cts {outFolder}/munged_ldsc_gwas/{gwas}_munged.sumstats.gz --ref-ld-chr /sc/arion/projects/ad-omics/ashvin/ldsc_annotations/1000G_EUR_Phase3_baseline/baseline. --out {outFolder}/enrichment_results/{gwas} --ref-ld-chr-cts {ldcts_prefix} --w-ld-chr /sc/arion/projects/ad-omics/ashvin/ldsc_annotations/weights_hm3_no_hla/weights."
