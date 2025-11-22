import sys 
import os.path
import glob
import pandas as pd 

gwas = config["GWAS"]

# this has to be fixed - ldsc requires 22 chroms
chromosome_numbers = list(range(1, 23))

# where the bed files are kept
infolder = config["inFolder"]

# list of file names
annot_prefix = config["prefix"]

#where to output the results
outFolder = config["outFolder"]
# revised_bed_file = ["revised_" + b for b in bed_file]

# flanking intervals
flank = config["flank"]

# file extension either .bed or .bed.gz
annot_ext = ".bed.gz"
annot_ext = config["annot_ext"]

ldcts_prefix = config["gwas_inputfile"]
gwas_df = pd.read_csv("/sc/arion/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary_GWAS.csv").set_index("dataset")
#gwas_df = pd.read_excel("/sc/arion/projects/ad-omics/data/references//GWAS/GWAS-QTL_data_dictionary.xlsx", sheet_name = 2).set_index("dataset")
genomeBuild = config['annotation_build']
category = ''.join(ldcts_prefix).split('.')[0]

print("test")
print(ldcts_prefix)
print(category)

# gwas_inputfile in config is the name of the file that contains lists of the different features

# currently config insists that each feature must be listed by name rather than as an external list of feature BEDs



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


# can this handle gzipped BED files? 
# currently writes out annot files to the current working directory rather than the outFolder
rule liftOver_bed:
    input:
        expand("{inFolder}/{{annotprefix}}{ext}", inFolder = infolder, annotprefix = annot_prefix, ext = annot_ext)
    output:
        expand("annotation_database/{c}/input_bed_files/{{annotprefix}}{ext}", c = category, annotprefix = annot_prefix, ext = annot_ext)
    shell:
        "ml liftover/09-Jul-2019;"
        "sh scripts/liftover_script.sh {infolder}/{wildcards.annotprefix}{annot_ext} annotation_database/{category}/input_bed_files/{wildcards.annotprefix}{annot_ext} {genomeBuild}"

# requires valr package to be installed on local R    
# here - can we add flanking of ranges?
# output isn't zipped so is huge for some reason    
rule create_annotation:
    input: 
        expand("annotation_database/{c}/input_bed_files/{{annotprefix}}{ext}", c = category, annotprefix = annot_prefix, ext = annot_ext)
    output: 
        expand("annotation_database/{c}/annotation_files/{{annotprefix}}.annot.gz", c = category, annotprefix = annot_prefix)
    params: 
        script1="scripts/create_annotation.R"
    shell: 
        "ml R;"
        "Rscript {params.script1} -b {input} -o {output} --flank {flank};"

# now the annotations get split
# appears to be run over and over rather than just producing 22 output files from one run
# why is this a separate script? couldn't this be done in the create_annotation step?
rule split_by_chr: 
    input: 
        expand("annotation_database/{c}/annotation_files/{{annotprefix}}.annot.gz", c = category, annotprefix = annot_prefix)
    output: 
        expand("annotation_database/{c}/annotation_files/{{annotprefix}}.a.{{chr_nums}}.annot.gz", c = category, chr_nums = chromosome_numbers, annotprefix = annot_prefix)
    params:
        c = category
    shell:
        "python scripts/split_by_chr.py --annotation {input} --output annotation_database/{params.c}/annotation_files/{wildcards.annotprefix} --chrom {wildcards.chr_nums}"

# looks like LD score gets calculated per chromosome
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

# does this run per chrom or all together?
rule new_ldcts_file:   
    input:
        expand("annotation_database/{c}/annotation_files/{annotprefix}.a.{chr_nums}.annot.gz", c = category, chr_nums = chromosome_numbers, annotprefix = annot_prefix),
        expand("annotation_database/{c}/annotation_files/{annotprefix}.a.{chr_nums}.l2.ldscore.gz", c = category, chr_nums = chromosome_numbers, outfolder = outFolder, annotprefix = annot_prefix), 
        expand("annotation_database/{c}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M", c = category, chr_nums = chromosome_numbers, annotprefix = annot_prefix),
        expand("annotation_database/{c}/annotation_files/{annotprefix}.a.{chr_nums}.l2.M_5_50", c = category, chr_nums = chromosome_numbers, outfolder = outFolder, annotprefix = annot_prefix)   
    output:
        ldcts = os.path.join(outFolder, ldcts_prefix)
        #expand("{ldctsprefix}", ldctsprefix = ldcts_prefix, outfolder = outFolder)
    params: 
        ldctsprefix = ldcts_prefix,
        outfolder = outFolder,
        script1 = "scripts/create_ldcts.R",
        annotprefix = annot_prefix,
        joined_bar=lambda w, input: ",".join(annot_prefix),  # ', input' was added
        c = category
    run:
        directory_path = os.path.join("annotation_database", category, "annotation_files")
        full_paths = [os.path.join(directory_path, fname + ".a.") for fname in annot_prefix]
        # Create DataFrame
        df = pd.DataFrame({
        'file_name': annot_prefix,
        'file_path': full_paths
        })
        # Write to TSV without column names (header=False)
        df.to_csv(output.ldcts, sep='\t', index=False, header=False)
        # "echo ${params.annotprefix}"
        #"ml R/4.0.3;"
        #"Rscript {params.script1} -i annotation_database/{params.c}/annotation_files/ -b {params.joined_bar} -o {params.ldctsprefix};"

## remove annotations with <100 SNPs
rule filter_ldcts:
    input: 
        ldcts = os.path.join(outFolder, ldcts_prefix)
    output:
        ldcts = os.path.join(outFolder, ldcts_prefix + ".filtered") 
    params:
        script = "scripts/filter_ldcts_annotations.py"
    shell:
        "python {params.script} {input.ldcts} {output.ldcts} 100 1e-10 8"

# this should only have to run once per GWAS
# why is this not just a separate python script?
rule format_sumstats:
    output:
        expand("{outfolder}/formatted_ldsc_gwas/{GWAS}.sumstats.gz", outfolder = outFolder, GWAS=gwas)
    params:
        out_folder = outFolder
    run:
        import numpy as np
        #gwas_dict = pd.read_excel("/sc/arion/projects/ad-omics/data/references//GWAS/GWAS-QTL_data_dictionary.xlsx", sheet_name = 2)
        gwas_dict = pd.read_csv("/sc/arion/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary_GWAS.csv")
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
            print(new_gwas_dict['full_effect'].str.contains("OR|or"))
            print(gwas_shortened.iloc[:,4], gwas_shortened.iloc[:,5])
            if (new_gwas_dict['full_effect'].str.contains("OR|or").values[0] == True):
                # gwas_shortened['Z'] = gwas_shortened.iloc[:,4] - 1  
                gwas_shortened['Z'] = np.log10(gwas_shortened.iloc[:,4]) / gwas_shortened.iloc[:,5]
            else:
                gwas_shortened['Z'] = gwas_shortened.iloc[:,4] / gwas_shortened.iloc[:,5]
            gwas_shortened = gwas_shortened.iloc[:, [0, 1, 2, 3, 6]]
            gwas_shortened.columns = ['ID', 'Allele1', 'Allele2', 'P.value', 'Z_score']
            gwas_shortened.dropna(inplace=True)
            print(gwas_shortened.head())
            return(gwas_shortened)
        print(gwas)
        for g in gwas:
            formatted_gwas = format_gwas(g)
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
        "python /hpc/packages/minerva-common/ldsc/1.0.0/ldsc/munge_sumstats.py "
        "--sumstats {wildcards.outfolder}/formatted_ldsc_gwas/{wildcards.GWAS}.sumstats.gz --merge-alleles  data/w_hm3.snplist "
        "--out {wildcards.outfolder}/munged_ldsc_gwas/{wildcards.GWAS}_munged "
        "--N {params.sample_size} --snp ID --a1 Allele1 --a2 Allele2 --p P.value --signed-sumstats Z_score,0"

# this is run once per GWAS / annotation pair, right?
rule run_ldsc:
    input:
        ldcts = os.path.join(outFolder, ldcts_prefix + ".filtered"),
        gwas = expand("{outfolder}/munged_ldsc_gwas/{{GWAS}}_munged.sumstats.gz", outfolder=outFolder, GWAS=gwas)
    output:
        expand("{outfolder}/enrichment_results/{{GWAS}}.cell_type_results.txt", outfolder=outFolder, GWAS=gwas)
    shell:
        "ml ldsc/1.0.0;"
        "python /hpc/packages/minerva-common/ldsc/1.0.0/ldsc/ldsc.py --h2-cts {outFolder}/munged_ldsc_gwas/{wildcards.GWAS}_munged.sumstats.gz "
        "--ref-ld-chr /sc/arion/projects/ad-omics/ashvin/ldsc_annotations/1000G_EUR_Phase3_baseline/baseline. --out {outFolder}/enrichment_results/{wildcards.GWAS} "
        "--ref-ld-chr-cts {input.ldcts} --w-ld-chr /sc/arion/projects/ad-omics/ashvin/ldsc_annotations/weights_hm3_no_hla/weights."
