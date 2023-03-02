suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(valr))
suppressPackageStartupMessages(library(optparse))
#!/usr/bin/env Rscript

thousand_genome_annotations = read_tsv('/sc/arion/projects/ad-omics/ashvin/ldsc/1000G_Phase3_annotation_full.tsv.gz')

create_annotation <- function(data_file) {
  print("decoy")
  bed_file = data_file
  df = thousand_genome_annotations
  colnames(bed_file) = c('chrom', 'start', 'end', 'category')
  cat = bed_file$category[1]
  print(head(cat))
  overlap = bed_intersect(df, bed_file)
  type_of_annotation = as.character(bed_file$category[1]) 
  df$temp <- 0 
  df$temp[thousand_genome_annotations$SNP %in% overlap$SNP.x] <- 1
  df$chrom = substr(df$chrom, 4,nchar(df$chrom))
  df$chrom = as.integer(df$chrom)
  df$start = as.integer(df$start)
  df = df[, c(1,2,3,4,5,7)]
  colnames(df) = c('CHR', 'BP', 'SNP', 'CM', 'base', 'ANNOT')
  df$ANNOT = as.integer(df$ANNOT)
  return(df)
}

option_list <- list(
   make_option(c('-b', '--bed'), help='list of bed files'),
   make_option(c('-o', '--out'), help='list of output files')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)
print(opt)

for (i in opt$bed) {
   bed_file_path <- i
   bed_file <- read.table(bed_file_path)
   output <- opt$out
   print(output)
   annotation_file <- create_annotation(bed_file)
   file_path = paste(output, '.annot', sep = '')
   print(file_path)
   write_tsv(annotation_file, output)
}
