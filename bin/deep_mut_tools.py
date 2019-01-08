#!/usr/bin/env python3
"""
Functions and Classes for processing Deep Mutational Scanning data
"""
import fileinput
import re
#import os
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
        """Generate a pandas series containing genotypes (as a list [X1Y, A2B, ...]) for the data"""
        geno = self.variant_data['variants'].str.replace('p.', '')
        geno = geno.str.split(',')
        if inc_wt:
            geno[geno.isna()] = 'WT'
        else:
            geno = geno.dropna()
        return geno

    def unique_variants(self):
        """Generate a list of the different variants in the dataset"""
        var = self.genotypes(inc_wt=False)
        var = list(set([i for x in var for i in x]))
        var.sort(key=lambda x: int(x[1:-1]))
        return var

    def write_ref_fasta(self, fasta_file):
        """Write gene reference sequence to a fasta file"""
        print(f">{self.meta_data['gene_name']}", file=fasta_file)
        print(split_lines(self.meta_data['ref_seq']), file=fasta_file)

    def write_variant_fasta(self, fasta_file):
        """Write a fasta file with all variant sequences"""
        gene = self.meta_data["gene_name"]
        uniprot_id = self.meta_data["uniprot_id"]

        for muts in self.genotypes():
            if muts == 'WT':
                print(f'>{gene}|{uniprot_id}|', file=fasta_file)
                print(split_lines(self.meta_data['ref_seq']), file=fasta_file, end='\n\n')

            else:
                seq = list(self.meta_data['ref_seq'])
                for var in i:
                    seq[int(var[1:-1]) - 1] = var[-1]

                print(f'>{gene}|{uniprot_id}|{','.join(muts)}', file=fasta_file)
                print(split_lines(''.join(seq)), file=fasta_file, end='\n\n')

def split_lines(seq, line_len=FA_LINE_LEN):
    """Split a string into lines of length line_len, inserting line breaks"""
    return '\n'.join([seq[i:i+line_len] for i in range(0, len(seq), line_len)])

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
