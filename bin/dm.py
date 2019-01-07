#!/usr/bin/env python3
"""
Template Script
"""
import sys
import os
import argparse
from contextlib import contextmanager
import deep_mut_tools as dm

class TaskSelecter:
    """Case selecter for different possible tasks to apply to DeepMut data"""
    def choose(self, task):
        """Choose the required task"""
        try:
            return getattr(self, task)
        except:
            raise ValueError(f"Task '{task}' not found")

    def ref_fasta(self, **kwargs):
        """Print reference sequence to out_file"""
        kwargs['deep_data'].write_ref_fasta(kwargs['out_file'])

    def variant_fasta(self, **kwargs):
        """Print all variant seqs to a fasta file"""
        kwargs['deep_data'].write_variant_fasta(kwargs['out_file'])

def main(args):
    """Main script"""
    deep_data = dm.read_deep_mut(args['dm_file'])

    if os.path.isdir(args['out']):
        raise ValueError('Output path is a directory')

    elif not os.path.exists(args['out']) or args['overwrite']:
        with open_file(args['out']) as out_file:
            selecter = TaskSelecter()
            selecter.choose(args['task'])(deep_data=deep_data, out_file=out_file, **args)

    else:
        raise FileExistsError(f"{args['out']} already exists")

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

    parser.add_argument('--out', '-o', default='', help="Output file, defaults to stdout")
    parser.add_argument('--overwrite', '-w', action='store_true',
                        help='Overwrite output file if it exists')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(vars(ARGS))
