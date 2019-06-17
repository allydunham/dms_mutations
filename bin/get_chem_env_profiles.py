#!/usr/bin/env python3
"""
Calculate chemical environment profiles for all positions in a PDB file(s) or PDB files
linked from dm files. If a single PDB file is given output will be to stdout unless a file
suffix is supplied. Otherwise output files are written in the directories of input
PDB files.
"""
import sys
import os
import argparse
import Bio.PDB
from Bio.SeqUtils import seq1
import chemical_environment as ce
import deep_mut_tools as dmt

# TODO don't assume model 0 always?
# TODO update to generate all chemical environment profiles selected

def main(args):
    """Main script"""

    # generate dictionary of profile functions
    profile_funcs = {}
    if args.k_nearest > 0:
        fun = lambda x, y, z: ce.k_nearest_residues_profile(x, y, z, args.k_nearest)
        profile_funcs[f'nearest_{args.k_nearest}'] = fun

    if args.angstroms > 0:
        fun = lambda x, y, z: ce.within_distance_profile(x, y, z, args.angstroms)
        profile_funcs[f'within_{args.angstroms}A'] = fun

    if not profile_funcs:
        raise ValueError(('No profile method selected, at least one of'
                          '--k_nearest/--angstroms must be used'))

    # construct list of PDB files to process
    if args.dm:
        pdb_files = []
        for dm_file in args.input:
            root = '/'.join(dm_file.split('/')[:-1])
            dm_header = dmt.read_deep_mut_header(dm_file)
            pdb_ids = [i.split(':') for i in dm_header['pdb_id']]

            for pdb in pdb_ids:
                pdb_files.append({'pdb_file': f'{root}/{pdb[0]}/{pdb[0]}.pdb',
                                  'chains': None if args.chains else [pdb[1]],
                                  'combine_chains': args.chains,
                                  'offset': pdb[2] if len(pdb) > 2 else 0})

    else:
        pdb_files = [{'pdb_file': f, 'chains': None,
                      'combine_chains': args.chains, 'offset': 0} for f in args.input]

    suffix = '.chem_env' if not args.suffix and len(pdb_files) > 1 else args.suffix

    # Calculate profiles for each PDB file
    pdb_parser = Bio.PDB.PDBParser()
    for pdb in pdb_files:
        structure = pdb_parser.get_structure(pdb['pdb_file'].split('/')[-1].split('.')[0],
                                             pdb['pdb_file'])
        residue_lists = make_residue_list(structure[0], chains=pdb['chains'],
                                          combine_chains=pdb['combine_chains'])

        pdb['residues'] = []
        pdb['profiles'] = {k: [] for k in profile_funcs}
        for residues in residue_lists:
            residues = ce.drop_hetero_atoms(residues)
            pdb['residues'].extend(residues)
            for key, func in profile_funcs.items():
                pdb['profiles'][key].extend(ce.get_profiles(residues,
                                                            func,
                                                            drop_hetero=False))

    # Write out profile files
    if len(pdb_files) == 1 and not args.dm and not args.suffix:
        write_profile_table(pdb_files[0]['residues'], pdb_files[0]['profiles'],
                            pdb_files[0]['offset'], file=sys.stdout)

    else:
        for pdb in pdb_files:
            with open(f"{os.path.splitext(pdb['pdb_file'])}{suffix}", 'w') as outfile:
                write_profile_table(pdb['residues'], pdb['profiles'],
                                    pdb['offset'], file=outfile)


def make_residue_list(model, chains=None, combine_chains=True):
    """
    Extract a list of residues lists from a Bio.PDB model.

    model: Bio.PDB model
    chain: list of chains to process, or None for all chains
    combine_chains: group chains together

    return: list of Bio.PDB residue lists
    """
    if combine_chains:
        return [[r for c in chains for r in model[c].get_residues()]]

    return [list(model[c].get_residues()) for c in chains]

# TODO map all profiles to strings before printing
def write_profile_table(residues, profiles, offset=0, file=sys.stdout):
    """
    Write a protein residue chemical environment profile table to file

    residues: list of Bio.PDB residues
    profiles: dict of chemical profiles (keys used as column names)
    offset: shift PDB sequence numbering (int)
    file: write to this file object

    returns: None
    """
    prof_keys = profiles.keys()
    print('pdb_id', 'chain', 'pos', 'aa', *prof_keys, sep='\t', file=file)
    for i, res in enumerate(residues):
        res_id = res.get_full_id()
        print(res_id[0], res_id[2], int(res_id[3][1]) + offset, seq1(res.get_resname()),
              *[','.join(map(str, profiles[k][i])) for k in prof_keys],
              sep='\t', file=file)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', metavar='I', nargs='+', help="Input file(s)")

    parser.add_argument('--dm', '-d', action='store_true', help="Assume input are dm files")

    parser.add_argument('--chains', '-c', action='store_true',
                        help="Process all chains together from PDB input")

    parser.add_argument('--suffix', '-s', default='',
                        help="Suffix to add to replace .pdb with in generated files")

    parser.add_argument('--k_nearest', '-k', default=0, type=int,
                        help='k for k-nearest amino acids profile')

    parser.add_argument('--angstroms', '-a', default=0, type=float,
                        help='Maximum distance (in A) to use in within distance profile')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
