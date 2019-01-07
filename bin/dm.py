#!/usr/bin/env python3
"""
Template Script
"""
import sys
import os
import argparse
from contextlib import contextmanager
import deep_mut_tools as dm

class DMTaskSelecter:
    """Case selecter for different possible tasks to apply to DeepMut data"""
    def __init__(self, deep_data, path):
        self.deep_data = deep_data
        self.path = path

    def choose(self, task):
        """Choose the required task"""
        try:
            return getattr(self, task)
        except:
            raise ValueError(f"Task '{task}' not found")

    def ref_fasta(self, **kwargs):
        """Print reference sequence to out_file"""
        with open_file(self.path) as out_file:
            self.deep_data.write_ref_fasta(out_file)

    def variant_fasta(self, **kwargs):
        """Print all variant seqs to a fasta file"""
        with open_file(self.path) as out_file:
            self.deep_data.write_variant_fasta(out_file)

    def sift4g(self, **kwargs):
        """Write the required files to fun sift4g on the data"""
        out_dir = self.path.rstrip('/') if self.path else '/'.join(self.path.split('/')[0:-1])
        gene_name = self.deep_data.meta_data['gene_name']

        fasta_path = f"{out_dir}/{gene_name}.fa"
        if not os.path.exists(fasta_path) or kwargs['overwrite']:
            with open_file(fasta_path) as fasta_file:
                self.deep_data.write_ref_fasta(fasta_file)

        # Export sorted list of variants
        variants = self.deep_data.variant_data.variants.str.split(',').dropna()
        variants = list(set([i.strip('p.') for x in variants for i in x]))
        variants.sort(key=lambda x: int(x[1:-1]))
        with open(f"{out_dir}/{gene_name}.subst", mode='w') as subst_file:
            print(*variants, sep='\n', file=subst_file)

def main(args):
    """Main script"""
    # Check if output file already exists and overwrite is not enabled
    # If it is a directory this is raises an error when trying to open it
    if not args['overwrite'] and os.path.isfile(args['out']):
        raise FileExistsError(f"File '{args['out']}' already exists and overwrite isn't enabled")

    deep_data = dm.read_deep_mut(args['dm_file'])
    selecter = DMTaskSelecter(deep_data, args['out'])
    selecter.choose(args['task'])(**args)

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

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('task', metavar='T', help="Task to perform on the deep mutagenesis file")
    parser.add_argument('dm_file', metavar='D', help="Deep mutagenesis file")

    parser.add_argument('--out', '-o', default='',
                        help="Output file or root directory (default: stdout or current directory)")
    parser.add_argument('--overwrite', '-w', action='store_true',
                        help='Overwrite output file if it exists')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(vars(ARGS))
