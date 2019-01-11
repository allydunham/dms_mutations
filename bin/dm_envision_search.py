#!/usr/bin/env python3
"""
Script to search envision
"""
import argparse
import deep_mut_tools as dm

def main(args):
    """Main script"""
    deep = dm.read_deep_mut(args.dm_path)

    # Extract genotypes in the format of Envision databases
    geno = deep.genotypes()
    geno = set([''.join([str(x) for x in var]) for sublist in geno for var in sublist])
    geno = [f'{deep.meta_data["uniprot_id"]}_{x}' for x in geno]

    with open(args.env_database, 'r') as env_file:
        # Fast seek to find protein first?

        # Slow see to find vars

        ### Use Grep instead? ###

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dm_path', metavar='D', help="dm file path")
    parser.add_argument('env_database', metavar='E', help="Envision database file path")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
