import pandas as pd
import numpy as np
import os
import time
import scipy.stats as stats
import logging
from pandas.api.types import is_integer_dtype
import gzip

def split_by_chr(annotation_file, name, output_path, chrom): 
    int_chrom = int(chrom)
    grouped_by_chrom = annotation_file[annotation_file['CHR'] == int_chrom]
    file_name = output_path + ".a." + chrom + ".annot.gz"
    print(chrom)
    print(grouped_by_chrom.head(5))
    grouped_by_chrom.to_csv(file_name, header=True, index=False, sep = '\t', compression='gzip')


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()    
    parser.add_argument('--annotation', required=True, help='Input annotation file')
    parser.add_argument('--output', required=True, help = 'Output annotation folder')
    parser.add_argument('--chrom', required=True, help = 'Chromosome number')
    args = parser.parse_args()
    print(args) 
    # read annotation file
    t0 = time.time()
    df_annotation = pd.read_csv(args.annotation, sep='\t')
    # print(df_annotation.head(5))
    df_annot = df_annotation
    #df_annot = df_annot[df_annot['CHR'] != 'X']
    #df_annot = df_annot[df_annot['CHR'] != 'Y']
    # split annotation file by chromosome
    df_annot['CHR'] = pd.to_numeric(df_annot['CHR'])
    # print(df_annot.head(5))
    split_by_chr(df_annot, args.annotation, args.output, args.chrom)

