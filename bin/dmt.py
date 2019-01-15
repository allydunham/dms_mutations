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
from ftplib import FTP
import gzip
import shutil
import pandas as pd
import evcouplings.utils as ev
import deep_mut_tools as dm
from nested_dicts import nested_merge
from smart_open import smart_open

ROTABASE_PATH = '/Users/ally/Projects/mutations/rotabase.txt'

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

        with smart_open(kwargs['path'], mode='w') as out_file:
            env.to_csv(out_file, index=False)

    def foldx(self, **kwargs):
        """Prepare a mutation list for FoldX analysis and fetch PDB file if it is not present"""
        path = kwargs['path']
        out_dir = path.rstrip('/') if path else '/'.join(kwargs['dm_file'].split('/')[0:-1])
        genotypes = self.deep_data.genotypes(wt=False, nonsense=False)

        for pdb in self.deep_data.meta_data['pdb_id']:
            # Download PDB if it doesn't exist
            pdb = pdb.split(':')
            pdb_id, pbd_chain = pdb[0], pdb[1]
            pdb_dir = f'{out_dir}/{pdb_id}'
            if not os.path.isdir(pdb_dir):
                os.mkdir(pdb_dir)

            with smart_open(f"{pdb_dir}/individual_list_{pdb_id}.txt", mode='w') as out_file:
                for geno in genotypes:
                    print(*[f"{i[0]}{pbd_chain}{i[1:]}" for i in geno],
                          sep=',', end=';\n', file=out_file)

            pdb_path = f"{pdb_dir}/{pdb_id}.pdb"
            pdb_gz_path = f"{pdb_path}.gz"
            pdb_id_low = pdb_id.lower()
            if not os.path.isfile(pdb_path):
                with open(pdb_gz_path, 'wb') as pdb_file, FTP('ftp.ebi.ac.uk') as ftp:
                    ftp.login('anonymous')
                    ftp.cwd(f'/pub/databases/pdb/data/structures/divided/pdb/{pdb_id_low[1:3]}/')
                    ftp.retrbinary(f'RETR pdb{pdb_id_low}.ent.gz', pdb_file.write)

                with gzip.open(pdb_gz_path, 'rb') as gz_file, open(pdb_path, 'wb') as pdb_file:
                    shutil.copyfileobj(gz_file, pdb_file)

                os.remove(pdb_gz_path)

            if not os.path.isfile(f'{pdb_dir}/rotabase.txt'):
                shutil.copy(kwargs['rotabase'], f'{pdb_dir}/rotabase.txt')

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

def main(args):
    """Main script"""
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

    parser.add_argument('--rotabase', '-r', help='Path to FoldX rotabase.txt file',
                        default=ROTABASE_PATH)

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(vars(ARGS))
