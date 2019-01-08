#!/usr/bin/env python3
"""
Tool to run various tasks on individual .dm files.
Currently implemented tasks:
- ref_fasta: generate a fasta file witht he ref seq
- variant_fasta: generate a fasta file with all variants
- sift4g: prepare the files to run sift4g
- evcouplings: prepare the files to run EVcouplings
"""
import sys
import os
import argparse
from contextlib import contextmanager
import pandas as pd
import evcouplings.utils as ev
import deep_mut_tools as dm
from nested_merge import nested_merge

class DMTaskSelecter:
    """Case selecter for different possible tasks to apply to DeepMut data"""
    def __init__(self, deep_data):
        self.deep_data = deep_data

    def choose(self, task):
        """Choose the required task"""
        try:
            return getattr(self, task)
        except:
            raise ValueError(f"Task '{task}' not found")

    def ref_fasta(self, **kwargs):
        """Print reference sequence to out_file"""
        with open_file(kwargs['path']) as out_file:
            self.deep_data.write_ref_fasta(out_file)

    def variant_fasta(self, **kwargs):
        """Print all variant seqs to a fasta file"""
        with open_file(kwargs['path']) as out_file:
            self.deep_data.write_variant_fasta(out_file)

    def sift4g(self, **kwargs):
        """Prepare the required files to fun sift4g on the data"""
        path = kwargs['path']
        out_dir = path.rstrip('/') if path else '/'.join(kwargs['dm_file'].split('/')[0:-1])
        gene_name = self.deep_data.meta_data['gene_name']
        fasta_path = f"{out_dir}/{gene_name}.fa"

        if kwargs['overwrite'] or not os.path.exists(fasta_path):
            self.ref_fasta(path=fasta_path, overwrite=False)

        # Export sorted list of variants
        variants = self.deep_data.unique_variants()
        with open(f"{out_dir}/{gene_name}.subst", mode='w') as subst_file:
            print(*variants, sep='\n', file=subst_file)

    def evcouplings(self, **kwargs):
        """Prepare the files required to run EVCouplings on the data"""
        path = kwargs['path']
        out_dir = path.rstrip('/') if path else '/'.join(kwargs['dm_file'].split('/')[0:-1])

        fasta_path = f"{out_dir}/{self.deep_data.meta_data['gene_name']}.fa"
        if kwargs['overwrite'] or not os.path.exists(fasta_path):
            self.ref_fasta(path=fasta_path, overwrite=False)

        csv_path = f"{out_dir}/ev_variants.csv"

        # Generate YAML config file
        if kwargs['ev_options']:
            overrides = ev.config.parse_config(kwargs['ev_options'])
        else:
            overrides = {}

        config = ev.config.read_config_file(kwargs['ev_default'])
        config = nested_merge(config,
                              {'global': {'prefix': f"{out_dir}/ev",
                                          'sequence_id': self.deep_data.meta_data['uniprot_id'],
                                          'sequence_file': fasta_path},
                               'mutate': {'mutation_dataset_file': csv_path}},
                              overrides)

        ev.config.write_config_file(f"{out_dir}/ev_config.txt", config)

        # Generate csv with variants
        variants = self.deep_data.variant_data[['variants', 'score']].copy()
        variants['variants'] = variants['variants'].str.replace('p.', '')
        variants.rename(index=str, columns={'variants':'mutants', 'score':'exp_score'})
        variants.to_csv(csv_path, sep=';', index=False)

    def envision(self, **kwargs):
        """Filter an envision database (arg env) to only include the desired variants"""
        if not kwargs['env']:
            raise ValueError("Must provide a precalculated envision database ('--env/-n') to "
                             "filter when performing the 'envision' task")

        env = pd.read_csv(kwargs['env'], comment='=').drop(columns='X1')
        env = env[env.Variant.isin(self.deep_data.unique_variants())]

        with open_file(kwargs['path'], 'w') as out_file:
            env.to_csv(out_file, index=False)


@contextmanager
def open_file(path, mode='w'):
    """Context manager managing a specified file or sys.stdout"""
    if not path:
        file_handle = sys.stdout
    else:
        file_handle = open(path, mode)
    yield file_handle
    if path:
        file_handle.close()

def main(args):
    """Main script"""
    # Check if output file already exists and overwrite is not enabled
    # If it is a directory this is raises an error when trying to open it
    deep_data = dm.read_deep_mut(args['dm_file'])
    selecter = DMTaskSelecter(deep_data)
    selecter.choose(args['task'])(**args)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('task', metavar='T', help="Task to perform on the deep mutagenesis file")
    parser.add_argument('dm_file', metavar='D', help="Deep mutagenesis file")

    parser.add_argument('--path', '-p', default='',
                        help="Output file/directory (default: stdout or input file directory)")
    parser.add_argument('--overwrite', '-o', action='store_true',
                        help='Overwrite generic files (e.g. ref fasta) during longer tasks '
                             'if they already exist')
    parser.add_argument('--ev_default', '-e', help='Base EVCouplings config file to modify',
                        default='/Users/ally/Projects/mutations/meta/base_evcouplings_config.txt')
    parser.add_argument('--ev_options', '-v', help='Additional EVCouplings config options')
    parser.add_argument('--env', '-n', help='Path to envision database file')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(vars(ARGS))
