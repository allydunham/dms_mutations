#!/usr/bin/env python3
"""
Functions and Classes for processing Deep Mutational Scanning data
"""
import fileinput
import re
import numpy as np
import pandas as pd

RE_FILE_VERSION = re.compile(r'#deep_mut_file_version')
RE_NUMERIC = re.compile(r'^-?[0-9]+(\.[0-9]+)?$')
RE_INT = re.compile(r'^-?[0-9]+$')

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
    t = read_deep_mut('data/standardised/araya_2012_hYAP65/variants.dm')
