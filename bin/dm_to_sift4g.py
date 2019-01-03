#!/usr/bin/env python3
"""
Convert a dm file to input for sift4g (a list of substitutions and a query fasta)
Output files are named gene_name.(fa/subst) based on the gene_name in the dm file, and will overwrite any existing files of that name.
"""
import argparse
import deep_mut_tools as dm

FA_LINE_LEN = 80

def main(args):
    """Main script"""
    deep_data = dm.read_deep_mut(args.dm_file)

    # Export fasta file
    with open(f"{deep_data.meta_data['gene_name']}.fa", 'w') as fasta_file:
        seq = deep_data.meta_data['ref_seq']
        print(f">{deep_data.meta_data['gene_name']}", file=fasta_file)
        for sub in [seq[i:i+FA_LINE_LEN] for i in range(0, len(seq), FA_LINE_LEN)]:
            print(sub, file=fasta_file)

    # Export list of variants
    variants = deep_data.variant_data.variants.str.split(',').dropna()
    variants = list(set([i.strip('p.') for x in variants for i in x]))
    variants.sort(key=lambda x: int(x[1:-1]))
    with open(f"{deep_data.meta_data['gene_name']}.subst", mode='w') as subst_file:
        print(*variants, sep='\n', file=subst_file)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dm_file', metavar='D', help="Input deep mutagenesis data (dm file)")

    parser.add_argument('--out', '-o', help='Output directory (default current dir)')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
