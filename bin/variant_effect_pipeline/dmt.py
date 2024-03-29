#!/usr/bin/env python3
"""
Tool to run various tasks on individual .dm files:
- ref_fasta: generate a fasta file witht he ref seq
- variant_fasta: generate a fasta file with all variants
- sift4g: prepare the files to run sift4g
- evcouplings: prepare the files to run EVcouplings
- polyphen2: prepare variants file for polyphen2 analysis
- foldx: prepare files for FoldX analysis of dm variants [NOT YET IMPLEMENTED]
"""
import os
import argparse
import gzip
import shutil
import collections
from ftplib import FTP

import pandas as pd
import evcouplings.utils as ev
from Bio.SeqUtils import seq1

import deep_mut_tools as dm
from nested_dicts import nested_merge
from smart_open import smart_open

ROTABASE_PATH = '/Users/ally/Projects/mutations/rotabase.txt'
FOLDX_JOB_SIZE = 1000

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
        self.deep_data.write_ref_fasta(kwargs['path'])

    def variant_fasta(self, **kwargs):
        """Print all variant seqs to a fasta file"""
        self.deep_data.write_variant_fasta(kwargs['path'])

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

        # Can pass the loaded config directly as a dict when using the class in other scripts
        try:
            config = ev.config.read_config_file(kwargs['ev_default'])
        except TypeError:
            config = kwargs['ev_default']

        config = nested_merge(config,
                              {'global': {'prefix': f"{out_dir}/ev",
                                          'sequence_id': self.deep_data.meta_data['uniprot_id'],
                                          'sequence_file': fasta_path},
                               'mutate': {'mutation_dataset_file': csv_path}},
                              overrides)

        ev.config.write_config_file(f"{out_dir}/ev_config.txt", config)

        # Generate csv with variants
        variants = self.deep_data.variant_data[['variants', 'score']].copy()
        variants.dropna(inplace=True)
        variants['variants'] = variants['variants'].str.replace('p.', '')
        variants = variants.rename({'variants':'mutant', 'score':'exp_score'}, axis='columns')
        variants.to_csv(csv_path, index=False)

    def envision(self, **kwargs):
        """Filter an envision database (arg env) to only include the desired variants"""
        if not kwargs['env']:
            raise ValueError("Must provide a precalculated envision database ('--env/-n') to "
                             "filter when performing the 'envision' task")

        env = pd.read_csv(kwargs['env'], comment='=').drop(columns='X1')
        env = env[env.Variant.isin(self.deep_data.unique_variants())]

        with smart_open(kwargs['path'], mode='w') as out_file:
            env.to_csv(out_file, index=False)

    def foldx(self, ftp=None, ftp_pdb_root='/pub/databases/pdb/data/structures/divided/pdb',
              **kwargs):
        """Prepare a mutation list for FoldX analysis and fetch PDB file if it is not present.
           Returns a dict giving the number of mutation list files generated for each PDB id"""
        path = kwargs['path']
        out_dir = path.rstrip('/') if path else '/'.join(kwargs['dm_file'].split('/')[0:-1])
        genotypes = self.deep_data.genotypes(wildtype=False, nonsense=False)

        generated_files = {}
        # Manage opening FTP connection, only if an opened FTP object has not been supplied
        # If a non EBI FTP connection is supplied it is expected the pdb_root
        # folder has the same format
        ftp_close = False
        ftp_pdb_root = '/pub/databases/pdb/data/structures/divided/pdb'
        if not isinstance(ftp, FTP):
            ftp = FTP('ftp.ebi.ac.uk')
            ftp.login('anonymous')
            ftp_close = True

        pdb_ids = self.deep_data.meta_data['pdb_id']
        for pdb in pdb_ids if isinstance(pdb_ids, list) else [pdb_ids]:
            # Download PDB if it doesn't exist
            pdb = pdb.split(':')
            pdb_id, chain = pdb[0], pdb[1]

            try:
                offset = int(pdb[2])
            except IndexError:
                offset = 0
            except ValueError as err:
                # If the given value is not a valid int
                raise err

            try:
                region = get_region(pdb[3])
            except IndexError:
                region = [0, float('inf')]
            except ValueError as err:
                # If the given value is not a valid region
                raise err

            pdb_dir = f'{out_dir}/{pdb_id}'
            if not os.path.isdir(pdb_dir):
                os.mkdir(pdb_dir)

            pdb_path = f"{pdb_dir}/{pdb_id}.pdb"
            pdb_gz_path = f"{pdb_path}.gz"
            pdb_id_low = pdb_id.lower()
            if not os.path.isfile(pdb_path):
                with open(pdb_gz_path, 'wb') as pdb_file:
                    ftp.cwd(f'{ftp_pdb_root}/{pdb_id_low[1:3]}/')
                    ftp.retrbinary(f'RETR pdb{pdb_id_low}.ent.gz', pdb_file.write)

                with gzip.open(pdb_gz_path, 'rb') as gz_file, open(pdb_path, 'wb') as pdb_file:
                    shutil.copyfileobj(gz_file, pdb_file)

                os.remove(pdb_gz_path)

            if not os.path.isfile(f'{pdb_dir}/rotabase.txt'):
                shutil.copy(kwargs['rotabase'], f'{pdb_dir}/rotabase.txt')

            subs = get_pdb_muts(pdb_path, single_letter=True)[chain]
            foldx_strs = [','.join(foldx_variants(g, subs, chain, offset, region)) for
                          g in genotypes]
            foldx_strs = [x for x in foldx_strs if x]

            if not kwargs['foldx_size']:
                with smart_open(f"{pdb_dir}/individual_list_{pdb_id}.txt", mode='w') as out_file:
                    print(*foldx_strs, sep=';\n', end=';\n', file=out_file)
                generated_files[pdb_id] = 1
            else:
                split_foldx_strs = [foldx_strs[i:i+kwargs['foldx_size']] for
                                    i in range(0, len(foldx_strs), kwargs['foldx_size'])]
                for i, var_set in enumerate(split_foldx_strs):
                    with smart_open(f"{pdb_dir}/individual_list_{pdb_id}_{i+1}.txt",
                                    mode='w') as out_file:
                        print(*var_set, sep=';\n', end=';\n', file=out_file)

                generated_files[pdb_id] = len(split_foldx_strs)

        if ftp_close:
            ftp.quit()

        return generated_files

    def polyphen2(self, **kwargs):
        """Generate protein variant file for Polyphen2 analysis"""
        variants = pd.DataFrame({'# Protein ID': self.deep_data.meta_data['uniprot_id'],
                                 'variant': self.deep_data.unique_variants()})

        variants['Position'] = variants.variant.str.slice(start=1, stop=-1)
        variants['AA1'] = variants.variant.str.slice(stop=1)
        variants['AA2'] = variants.variant.str.slice(start=-1)

        with smart_open(kwargs['path'], mode='w') as out_file:
            variants.to_csv(out_file, index=False, sep='\t',
                            columns=['# Protein ID', 'Position', 'AA1', 'AA2'])

def get_pdb_muts(path_or_file, single_letter=False):
    """Return a nested dictionary structure of sequence variants in a pdb file
       in the format x[chain][pdb residue num] = (ref, alt)"""
    with smart_open(path_or_file, 'r') as pdb_file:
        func = seq1 if single_letter else lambda x: x
        muts = collections.defaultdict(dict)
        for line in pdb_file:
            # Currently only deals with these conflict types
            # but can be extended to others if needed
            if line[:6] == 'SEQADV' and line[49:].lower().strip() in ('engineered',
                                                                      'variant',
                                                                      'microheterogeneity'):
                val = {'pdb_name': line[12:15].strip(),
                       'pdb_pos': int(line[18:22].strip()),
                       'chain': line[16:17].strip(),
                       'db_name': line[39:42].strip(),
                       'db_pos': int(line[43:48].strip()),
                       'conflict': line[49:].lower().strip()}

                muts[val['chain']][val['pdb_pos']] = (func(val['db_name']), func(val['pdb_name']))

    return muts

def foldx_variants(genotypes, sub, chain='A', offset=0, region=None):
    """Convert a list of genotypes into a list of FoldX individual_variants.txt entries.
       Returns an empty list if no variants are suitable (e.g. out of region or mutate to same)"""
    if region is None:
        region = [0, float('inf')]

    foldx_strs = []
    for geno in genotypes:
        ref, pos, mut = geno[0], int(geno[1:-1])-offset, geno[-1]
        if region[0] <= pos <= region[1] and not ref == mut:
            foldx_strs.append(f'{sub[pos][1] if pos in sub else ref}{chain}{pos}{mut}')

    return foldx_strs

def get_region(reg_str):
    """Transform a string of the form X-Y into a region"""
    spl = reg_str.split('-')
    if spl[0] == '' and spl[1] == '':
        # Fully open region
        reg = [0, float('inf')]
    elif spl[0] == '':
        # Open beginning
        reg = [0, int(spl[1])]
    elif spl[1] == '':
        # open end
        reg = [int(spl[0]), float('inf')]
    else:
        reg = [int(spl[0]), int(spl[1])]

    return reg

def main(args):
    """Main script"""
    deep_data = dm.read_deep_mut(args['dm_file'])
    selecter = DMTaskSelecter(deep_data)
    _ = selecter.choose(args['task'])(**args)

def non_neg_int(value):
    """Convert to non-negative integer or raise an error"""
    value = int(value)
    if value < 0:
        raise argparse.ArgumentTypeError(f'foldx_size ({value}) must be a non-negative integer')
    return value

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

    parser.add_argument('--foldx_size', '-f', type=non_neg_int, default=FOLDX_JOB_SIZE,
                        help='Number of variants to put in each FoldX Job (0 means all)')
    parser.add_argument('--rotabase', '-r', help='Path to FoldX rotabase.txt file',
                        default=ROTABASE_PATH)

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(vars(ARGS))
