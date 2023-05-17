library(readr)
library(stringr)
library(optparse)

#print("hello!")
create_ldcts_file <- function(bed_file_list, input_directory) {
   bed_file_list <- str_split(bed_file_list, "\\,", simplify=T)[,]
   #print(bed_file_list)
   data <- data.frame(bed_file_list)
   data$path <- c(paste(input_directory, bed_file_list, '.a.', sep = ''))
   return(data)
}


option_list <- list(
  make_option(c('-i', '--input'), help='input file directory'),
  make_option(c('-b', '--bed'), help='list of bed files'),
  make_option(c('-o', '--out'), help='name of ldcts file')
)


option.parser <- OptionParser(option_list=option_list)
arguments <- parse_args (option.parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args


# parser <-OptionParser(option_list=option_list)

#args <- commandArgs(trailingOnly = TRUE)


print(args)

myfilelist <- strsplit(opt$b, " ")
print(myfilelist)
ldcts_file <- create_ldcts_file(opt$bed, opt$input)

head(ldcts_file)
write.table(ldcts_file, opt$out, sep = '\t', row.names =FALSE, col.names = FALSE, quote = FALSE)
