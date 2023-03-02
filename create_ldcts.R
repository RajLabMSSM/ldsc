library(readr)
library(stringr)
library(optparse)

#print("hello!")
create_ldcts_file <- function(bed_file_list) {
   bed_file_list <- str_split(bed_file_list, "\\.", simplify=T)[,1]
   data <- data.frame(bed_file_list)
   data$path <- c(paste('bed_files/annotation_files/', bed_file_list, '.bed.a.', sep = ''))
   return(data)
}

option_list <- list(
  make_option(c('-b', '--bed'), help='list of bed files'),
  make_option(c('-o', '--out'), help='name of ldcts file')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)
#print(opt)
#print(opt$bed)

ldcts_file <- create_ldcts_file(opt$bed)
head(ldcts_file)
write.table(ldcts_file, opt$out, sep = '\t', row.names =FALSE, col.names = FALSE)
