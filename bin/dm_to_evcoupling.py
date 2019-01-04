#!/usr/bin/env python3
"""
Prepare data from a dm file for EVcoupling processing, creating the required YAML config file
"""
import argparse
import evcouplings.utils as ev
import deep_mut_tools as dm

def main(args):
    """Main script"""
    out_dir = args.out.rstrip('/') if args.out else '/'.join(args.out.split('/')[0:-1])
    deep_data = dm.read_deep_mut(args.deep_data)

    # Generate YAML config file
    config = ev.config.read_config_file(args.base)
    config['global']['prefix'] = f"{out_dir}/evcoup"
    config['global']['sequence_id'] = deep_data.meta_data['uniprot_id']

    fasta_path = f"{out_dir}/{deep_data.meta_data['gene_name']}.fa"
    deep_data.write_ref_fasta(path=fasta_path, overwrite=False)
    config['global']['sequence_file'] = fasta_path

    ev.config.write_config_file(f"{out_dir}/evcoup_config.txt", config)

    # Generate csv with variants
    deep_data.variant_data['variants'] = deep_data.variant_data['variants'].str.replace('p.', '')
    deep_data.variant_data = deep_data.variant_data[['variants', 'score']]
    deep_data.variant_data.rename(index=str, columns={'variants':'mutants', 'score':'exp_score'})
    deep_data.variant_data.to_csv(f"{out_dir}/ev_variants.csv", sep=';', index=False)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('deep_data', metavar='D', help="Input dm file")

    parser.add_argument('--out', '-o', default='',
                        help="Output directory, defaults to same as input file. It is used for config files and evcoupling output.")

    parser.add_argument('--base', '-b', help='Base config file to modify',
                        default='/Users/ally/Projects/mutations/meta/base_evcouplings_config.txt')


    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
