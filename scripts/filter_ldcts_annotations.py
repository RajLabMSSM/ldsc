import pandas as pd
import gzip
import os
import sys
from pathlib import Path

def check_annotation(annot_prefix, chromosomes=range(1, 23), min_snps=100, min_var=1e-10):
    """
    Check if an annotation is valid for LDSC across all chromosomes.
    
    Parameters:
    -----------
    annot_prefix : str
        Path prefix for annotation files (without chromosome number)
    chromosomes : range or list
        Chromosomes to check
    min_snps : int
        Minimum number of annotated SNPs required
    min_var : float
        Minimum variance required
    
    Returns:
    --------
    tuple : (is_valid, reason, stats_dict)
    """
    total_snps = 0
    total_annotated = 0
    all_vars = []
    missing_files = []
    
    for chrom in chromosomes:
        annot_file = f"{annot_prefix}{chrom}.annot.gz"
        
        if not os.path.exists(annot_file):
            missing_files.append(chrom)
            continue
            
        try:
            df = pd.read_csv(annot_file, sep='\t', compression='gzip')
            annot_col = df.iloc[:, -1]
            
            total_snps += len(df)
            total_annotated += (annot_col > 0).sum()
            all_vars.append(annot_col.var())
            
        except Exception as e:
            return False, f"Error reading chr{chrom}: {str(e)}", {}
    
    if missing_files:
        return False, f"Missing files for chromosomes: {missing_files}", {}
    
    if total_annotated < min_snps:
        return False, f"Too few annotated SNPs: {total_annotated} < {min_snps}", {
            'total_snps': total_snps,
            'annotated_snps': total_annotated,
            'mean_variance': sum(all_vars) / len(all_vars) if all_vars else 0
        }
    
    mean_var = sum(all_vars) / len(all_vars) if all_vars else 0
    if mean_var < min_var:
        return False, f"Variance too low: {mean_var} < {min_var}", {
            'total_snps': total_snps,
            'annotated_snps': total_annotated,
            'mean_variance': mean_var
        }
    
    return True, "Valid", {
        'total_snps': total_snps,
        'annotated_snps': total_annotated,
        'mean_variance': mean_var
    }

def filter_ldcts_file(input_ldcts, output_ldcts, min_snps=100, min_var=1e-10, verbose=True):
    """
    Filter an LDCTS file to remove invalid annotations.
    
    Parameters:
    -----------
    input_ldcts : str
        Path to input .ldcts file
    output_ldcts : str
        Path to output filtered .ldcts file
    min_snps : int
        Minimum number of annotated SNPs required
    min_var : float
        Minimum variance required
    verbose : bool
        Print progress information
    """
    valid_lines = []
    invalid_annotations = []
    
    print(f"Reading {input_ldcts}...")
    with open(input_ldcts, 'r') as f:
        lines = f.readlines()
    
    print(f"Checking {len(lines)} annotations...\n")
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            valid_lines.append(line)
            continue
        
        # Parse the line
        parts = line.split('\t')
        if len(parts) < 2:
            print(f"Warning: Skipping malformed line: {line}")
            continue
        
        annot_name = parts[0]
        annot_prefix = parts[1]
        
        if verbose:
            print(f"[{i+1}/{len(lines)}] Checking {annot_name}...", end=' ')
        
        # Check if annotation is valid
        is_valid, reason, stats = check_annotation(annot_prefix, min_snps=min_snps, min_var=min_var)
        
        if is_valid:
            valid_lines.append(line)
            if verbose:
                print(f"✓ VALID (SNPs: {stats['annotated_snps']}, Var: {stats['mean_variance']:.2e})")
        else:
            invalid_annotations.append((annot_name, reason, stats))
            if verbose:
                print(f"✗ INVALID - {reason}")
    
    # Write filtered file
    print(f"\nWriting filtered annotations to {output_ldcts}...")
    with open(output_ldcts, 'w') as f:
        for line in valid_lines:
            f.write(line + '\n')
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"SUMMARY:")
    print(f"{'='*60}")
    print(f"Total annotations: {len(lines)}")
    print(f"Valid annotations: {len(valid_lines)}")
    print(f"Invalid annotations: {len(invalid_annotations)}")
    print(f"Output written to: {output_ldcts}")
    
    if invalid_annotations:
        print(f"\n{'='*60}")
        print(f"INVALID ANNOTATIONS:")
        print(f"{'='*60}")
        for name, reason, stats in invalid_annotations:
            print(f"\n{name}:")
            print(f"  Reason: {reason}")
            if stats:
                print(f"  Stats: {stats}")
    
    return len(valid_lines), len(invalid_annotations)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python filter_ldcts_annotations.py <input.ldcts> <output.ldcts> [min_snps] [min_var]")
        print("\nExample:")
        print("  python filter_ldcts_annotations.py ENCODE.ldcts ENCODE_filtered.ldcts 100 1e-10")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    min_snps = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    min_var = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-10
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found!")
        sys.exit(1)
    
    print(f"Filtering parameters:")
    print(f"  Minimum annotated SNPs: {min_snps}")
    print(f"  Minimum variance: {min_var}")
    print()
    
    filter_ldcts_file(input_file, output_file, min_snps=min_snps, min_var=min_var)
