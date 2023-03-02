suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


#!/usr/bin/env Rscript

gwas_dict <- readxl::read_excel("/sc/arion/projects/ad-omics/data/references//GWAS/GWAS-QTL_data_dictionary.xlsx", sheet = 3)

format_gwas <- function(gwas) {
  new_gwas_dict <- gwas_dict[gwas_dict$dataset == gwas,]
  gwas_df <- read_tsv(new_gwas_dict$full_processed_path)
  column_names <- c(new_gwas_dict$full_snp, 
                    new_gwas_dict$full_A1, 
                    new_gwas_dict$full_A2, 
                    new_gwas_dict$full_p,
                    new_gwas_dict$full_effect,
                    new_gwas_dict$full_se)
#  print(column_names)
#  head(new_gwas_dict)
  gwas_shortened <- gwas_df[, column_names]
  if (grepl(new_gwas_dict$full_effect, "OR", fixed=TRUE)) {
    gwas_shortened$Z <- log(gwas_shortened[,5]) / gwas_shortened[,6]
  } else { 
    gwas_shortened$Z <- gwas_shortened[,5] / gwas_shortened[,6]
  }
  gwas_shortened <- gwas_shortened[, c(1,2,3,4,7)]
  colnames(gwas_shortened) = c('ID', 'Allele1', 'Allele2', 'P.value', 'Z_score')
  print(head(gwas_shortened))
  return(gwas_shortened)
}

option_list <- list(
  make_option(c('-g', '--gwas'), help='list of GWAS files'),
  make_option(c('-o', '--out'), help='list of output files')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)
print(opt)

formatted_gwas <- format_gwas(opt$gwas)
print("hello!")
setwd('/sc/arion/projects/ad-omics/ashvin/ldsc/formatted_ldsc_gwas')
write_tsv(formatted_gwas, opt$out)
