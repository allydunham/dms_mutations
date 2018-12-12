#!/usr/bin/env python3
"""
Functions and Classes for processing Deep Mutational Scanning data
"""
import fileinput
import numpy as np
import pandas as pd

class DeepMut:
    """
    Deep Mutational Scanning data with associated metadata
    """
    required_metadata = {'gene_name':None, 'domain':None, 'species':None, 'ref_seq':None,
                         'transform':'None', 'uniprot_id':None, 'authour':None, 'year':None}

    def __init__(self, variant_data, **kwargs):
        # Simple tests for variant data format
        if not isinstance(variant_data, 'DataFrame'):
            raise TypeError('variant_data must of type DataFrame (see Pandas)')
        elif not all([x in variant_data.columns for x in ('variants', 'score', 'raw_score')]):
            raise ValueError('variant_data must have columns: "variants", "score", "raw_score"')

        self.variant_data = variant_data

        # Set meta data and process any missing required fields
        self.meta_data = kwargs
        for key, value in self.required_metadata.items():
            if not key in self.meta_data.keys():
                self.meta_data[key] = value

def read_deep_mut(path):
    """Import deep mutational scanning data from '.dm' file and return a DeepMut object"""
    variant_data = pd.read_csv(path, sep='\t', comment='#')

    meta = {}
    with fileinput.input(path) as file_obj:
        for line in file_obj:
            if line[0] == '#':
                pass
            elif line[0] == '?':
                break
            else:
                raise ValueError('Reached the end of metadata lines (marked by #) before'
                                 'encountering the variant data table header (marked ?)')


if __name__ == "__main__":
    pass
