#!/usr/bin/env python3
"""
Functions and Classes for processing Deep Mutational Scanning data
"""
import fileinput
import re
#import numpy as np
import pandas as pd

RE_FILE_VERSION = re.compile(r'#deep_mut_file_version')
RE_NUMERIC = re.compile(r'^-?[0-9]+(\.[0-9]+)?$')
RE_INT = re.compile(r'^-?[0-9]+$')
FA_LINE_LEN = 80

class DeepMut:
    """
    Deep Mutational Scanning data with associated metadata
    """
    required_metadata = {'gene_name':None, 'domain':None, 'species':None, 'ref_seq':None,
                         'transform':'None', 'uniprot_id':None, 'authour':None, 'year':None}

    def __init__(self, variant_data, **kwargs):
        # Simple tests for variant data format
        if not isinstance(variant_data, pd.core.frame.DataFrame):
            raise TypeError('variant_data must of type DataFrame (see Pandas)')
        elif not all([x in variant_data.columns for x in ('variants', 'score', 'raw_score')]):
            raise ValueError('variant_data must have columns: "variants", "score", "raw_score"')

        self.variant_data = variant_data

        # Set meta data and process any missing required fields
        self.meta_data = kwargs
        for key, value in self.required_metadata.items():
            if not key in self.meta_data.keys():
                self.meta_data[key] = value

    def print_head(self, num=5):
        """Print meta_data and the head of variant_data in a readable format"""
        for key, value in self.meta_data.items():
            print(f'{key}:{value}')

        print(self.variant_data.head(num))

    def genotypes(self, inc_wt=True):
        """Generate a list of genotypes for the data"""
        geno = []
        for i in self.variant_data.variants:
            if not pd.isna(i):
                i = [x.strip('p.') for x in i.split(',')]
                geno.append([(x[0], int(x[1:-1]), x[-1]) for x in i])
            elif inc_wt:
                # Assume empty variants field means wt
                geno.append('WT')
        return geno

    def write_ref_fasta(self, path=""):
        """Write gene reference sequence to a fasta file"""
        if not path:
            path = f"{self.meta_data['gene_name']}.fa"

        with open(path, 'w') as fasta_file:
            seq = self.meta_data['ref_seq']
            print(f">{self.meta_data['gene_name']}", file=fasta_file)
            for sub in [seq[i:i+FA_LINE_LEN] for i in range(0, len(seq), FA_LINE_LEN)]:
                print(sub, file=fasta_file)




def read_deep_mut(path):
    """Import deep mutational scanning data from '.dm' file and return a DeepMut object"""
    # Import variants
    variant_data = pd.read_csv(path, sep='\t', comment='#')
    variant_data.columns = [x.strip('?') for x in variant_data.columns]

    # Import meta data
    meta = {}
    with fileinput.input(path) as file_obj:
        for line in file_obj:
            if RE_FILE_VERSION.match(line):
                continue

            if line[0] == '#' and not RE_FILE_VERSION.match(line):
                key, value = line.strip('#\n\t ').split(':', 1)
                if key == 'ref_seq':
                    seq = []
                    for _ in range(int(value)):
                        seq.append(file_obj.readline().strip('#+\n\t '))
                    meta['ref_seq'] = ''.join(seq)

                elif RE_INT.match(value):
                    meta[key] = int(value)

                elif RE_NUMERIC.match(value):
                    meta[key] = float(value)

                else:
                    meta[key] = value

            elif line[0] == '?':
                break

            else:
                raise ValueError('Reached the end of metadata lines (marked by #) before '
                                 'encountering the variant data table header (marked ?)')

    return DeepMut(variant_data, **meta)


if __name__ == "__main__":
    TEST = read_deep_mut('data/standardised/araya_2012_hYAP65/variants.dm')
