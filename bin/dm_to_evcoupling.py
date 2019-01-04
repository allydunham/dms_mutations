#!/usr/bin/env python3
"""
Prepare data from a dm file for EVcoupling processing, creating the required YAML config file
"""
import argparse
import evcouplings as ev
import deep_mut_tools as dm

def main(args):
    """Main script"""
    out_dir = args.out.rstrip('/')
    deep_mut = dm.read_deep_mut(args.deep_data)

    # Generate YAML config file

    # Generate csv with variants
    deep_mut.variant_data['variants'] = deep_mut.variant_data['variants'].str.replace('p.', '')
    deep_mut.variant_data = deep_mut.variant_data[['variants', 'score']]
    deep_mut.variant_data.rename(index=str, columns={'variants':'mutants', 'score':'exp_score'})
    deep_mut.variant_data.to_csv(f"{out_dir}/ev_variants.csv", sep=';', index=False)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('deep_data', metavar='D', help="Input dm file")

    parser.add_argument('--out', '-o', default='.', help="Output directory")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
