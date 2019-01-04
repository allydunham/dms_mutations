#!/usr/bin/env python3
"""
Convert a dm file to input for sift4g (a list of substitutions and a query fasta)
Output files are named gene_name.(fa/subst) based on the gene_name in the dm file, and will overwrite any existing files of that name.
"""
import argparse
import os
import deep_mut_tools as dm

def main(args):
    """Main script"""
    deep_data = dm.read_deep_mut(args.dm_file)
    out_dir = args.out.rstrip('/') + '/'

    # Export fasta file if it does not already exist (assumes a correctly named .fa file is right)
    fasta_path = f"{out_dir}{deep_data.meta_data['gene_name']}.fa"
    deep_data.write_ref_fasta(path=fasta_path, overwrite=False)

    # Export list of variants
    variants = deep_data.variant_data.variants.str.split(',').dropna()
    variants = list(set([i.strip('p.') for x in variants for i in x]))
    variants.sort(key=lambda x: int(x[1:-1]))
    with open(f"{out_dir}{deep_data.meta_data['gene_name']}.subst", mode='w') as subst_file:
        print(*variants, sep='\n', file=subst_file)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dm_file', metavar='D', help="Input deep mutagenesis data (dm file)")

    parser.add_argument('--out', '-o', default='.', help='Output directory (default current dir)')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
