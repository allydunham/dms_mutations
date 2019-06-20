#!/usr/bin/env python3
"""
Functions and Classes for processing Deep Mutational Scanning data
"""
import fileinput
import re
import copy
import pandas as pd
from smart_open import smart_open

# TODO enforce expected types of fields

DM_FILE_VERSION = '2.1.0'
RE_NUMERIC = re.compile(r'^-?[0-9]+(\.[0-9]+)?$')
RE_INT = re.compile(r'^-?[0-9]+$')
FA_LINE_LEN = 80

class DeepMut:
    """
    Deep Mutational Scanning data with associated metadata
    """
    required_metadata = {'gene_name':None, 'domain':None, 'species':None, 'aa_seq':None,
                         'transform':'None', 'uniprot_id':None, 'authour':None, 'year':None,
                         'pdb_id':None}

    # Used for grouping when writing files, to make them more readable
    key_groups = {'gene_keys': ['gene_name', 'domain', 'species', 'alt_name'],
                  'study_keys': ['authour', 'year', 'title', 'pmid', 'url', 'doi']}

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

    def genotypes(self, wildtype=True, nonsense=True):
        """Generate a pandas series containing genotypes (as a list [X1Y, A2B, ...]) for the data"""
        geno = self.variant_data['variants'].str.replace('p.', '')
        if wildtype:
            geno[geno.isna()] = 'WT'
        else:
            geno = geno.dropna()

        if not nonsense:
            geno = geno.dropna()
            geno = geno[~geno.str.contains(r'\*')]

        geno = geno.str.split(',')
        return geno

    def unique_variants(self, nonsense=True):
        """Generate a list of the different variants in the dataset"""
        var = self.genotypes(wildtype=False, nonsense=nonsense)
        var = list(set([i for x in var for i in x]))
        var.sort(key=lambda x: int(x[1:-1]))
        return var

    def write_ref_fasta(self, path_or_file='', mode='w', seq='aa_seq'):
        """Write gene reference sequence to a fasta file"""
        with smart_open(path_or_file, mode=mode) as fasta_file:
            print(f">{self.meta_data['gene_name']}", file=fasta_file)
            print(*split_lines(self.meta_data[seq]), sep='\n', file=fasta_file)

    def write_variant_fasta(self, path_or_file='', mode='w'):
        """Write a fasta file with all variant sequences"""
        gene = self.meta_data["gene_name"]
        uniprot_id = self.meta_data["uniprot_id"]

        with smart_open(path_or_file, mode=mode)as fasta_file:
            for muts in self.genotypes():
                if muts == 'WT':
                    print(f'>{gene}|{uniprot_id}|', file=fasta_file)
                    print(*split_lines(self.meta_data['aa_seq']), sep='\n',
                          file=fasta_file, end='\n\n')

                else:
                    seq = list(self.meta_data['aa_seq'])
                    for var in muts:
                        seq[int(var[1:-1]) - 1] = var[-1]

                    mut_str = ','.join(muts)
                    print(f">{gene}|{uniprot_id}|{mut_str}", file=fasta_file)
                    print(*split_lines(''.join(seq)), sep='\n', file=fasta_file, end='\n\n')

    def write_dm(self, path_or_file='', mode='w'):
        """Write data to dm file"""
        # Determine ordering of keys for writing meta data
        self.key_groups['accession_keys'] = [x for x in self.meta_data if '_id' in x]
        self.key_groups['sequence_keys'] = [x for x in self.meta_data if '_seq' in x]

        self.key_groups['misc_keys'] = [x for x in self.meta_data if x not in
                                        [x for _, v in self.key_groups.items() for x in v] and
                                        not x == 'deep_mut_file_version']

        ordered_keys = self.key_groups['gene_keys'] +\
                       self.key_groups['accession_keys'] +\
                       self.key_groups['study_keys'] +\
                       self.key_groups['misc_keys'] +\
                        self.key_groups['sequence_keys']



        with smart_open(path_or_file, mode=mode) as dm_file:
            print(f'#deep_mut_file_version:{DM_FILE_VERSION}', file=dm_file)

            for key in ordered_keys:
                if self.meta_data[key] is None or not key in self.meta_data:
                    print(f"#{key}:NA", file=dm_file)
                elif isinstance(self.meta_data[key], (list, tuple, set)):
                    lines = ['#*' + x for x in self.meta_data[key]]
                    print(f"#*{key}:{len(lines)}", file=dm_file)
                    print(*lines, sep='\n', file=dm_file)
                elif '_seq' in key:
                    lines = ['#+' + x for x in split_lines(self.meta_data[key])]
                    print(f"#+{key}:{len(lines)}", file=dm_file)
                    print(*lines, sep='\n', file=dm_file)
                else:
                    print(f"#{key}:{self.meta_data[key]}", file=dm_file)

            print('?', '\t'.join(self.variant_data.columns.values), sep='', file=dm_file)
            self.variant_data.to_csv(dm_file, header=False, index=False, sep='\t')


def split_lines(seq, line_len=FA_LINE_LEN):
    """Split a string into lines of length line_len, inserting line breaks"""
    return [seq[i:i+line_len] for i in range(0, len(seq), line_len)]

def read_deep_mut(path):
    """Import deep mutational scanning data from '.dm' file and return a DeepMut object"""
    # Import variants
    variant_data = pd.read_csv(path, sep='\t', comment='#')
    variant_data.columns = [x.strip('?') for x in variant_data.columns]

    # Import meta data
    meta = read_deep_mut_header(path)
    return DeepMut(variant_data, **meta)

def read_deep_mut_header(path):
    """Import header of a .dm file as a dictionary"""
    meta = copy.copy(DeepMut.required_metadata)
    with fileinput.input(path) as file_obj:
        for line in file_obj:
            if line[0] == '#':
                key, value = line.strip('#*+\n\t ').split(':', 1)
                if line[1] in ('*', '+'):
                    # read lists and extended entries
                    meta[key] = []
                    for _ in range(int(value)):
                        meta[key].append(file_obj.readline().strip('#+*\n\t '))

                    if line[1] == '+':
                        meta[key] = ''.join(meta[key])

                else:
                    # read normal entries
                    meta[key] = value

                meta[key] = weak_numeric_conversion(meta[key])

            elif line[0] == '?':
                break

            else:
                raise ValueError('Reached the end of metadata lines (marked by #) before '
                                 'encountering the variant data table header (marked ?)')
    return meta

def weak_numeric_conversion(var):
    """Convert x into ints or floats if appropriate or return x
       If x is a list, tuple or set a list is returned
       if a string an int or float
       otherwise x is returned unchanged"""
    if isinstance(var, (str, bytes)):
        if RE_INT.match(var):
            var = int(var)

        elif RE_NUMERIC.match(var):
            var = float(var)

    elif isinstance(var, (list, tuple, set)):
        if all(RE_INT.match(i) for i in var):
            var = [int(i) for i in var]

        elif all(RE_NUMERIC.match(i) for i in var):
            var = [float(i) for i in var]

    return var



if __name__ == "__main__":
    TEST = read_deep_mut('/Users/ally/Projects/mutations/data/standardised/'
                         'araya_2012_hYAP65/P46937_YAP1.dm')
