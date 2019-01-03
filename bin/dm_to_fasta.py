#!/usr/bin/env python3
"""
Script to convert a dm file into a fasta file with one sequence per variant
"""
import argparse
import deep_mut_tools as dm

def main(args):
    """Main script"""
    deep = dm.read_deep_mut(args.dm_path)
    geno = deep.genotypes()
    for i in geno:
        if i == 'WT':
            print(f'>{deep.meta_data["gene_name"]}|{deep.meta_data["uniprot_id"]}|')
            print(deep.meta_data['ref_seq'])

        else:
            # Print header
            muts = ','.join([''.join([str(y) for y in x]) for x in i])
            print(f'>{deep.meta_data["gene_name"]}|{deep.meta_data["uniprot_id"]}|{muts}')

            # Print seq
            seq = list(deep.meta_data['ref_seq'])
            for var in i:
                seq[var[1] - 1] = var[2]
            print(''.join(seq), '\n')

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dm_path', metavar='D', help="Path to input dm file")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
