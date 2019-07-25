#!/usr/bin/env python3
"""
Combine FoldX output files into a single table. Done by crawling output directories of
structure Uniprot ID/PDB ID/FoldX Output. Outputs to stdout.
"""

import sys
import os
import argparse
import pandas as pd

def main(args):
    """
    Main
    """
    write_header = True
    for top_dir in args.dirs:
        for path, _, filelist in os.walk(top_dir):
            if not any('.fxout' in i for i in filelist):
                continue

            uniprot_id, pdb_chain = path.split('/')[-2:]
            pdb_id, chain = pdb_chain.split('_')

            fx_output_file = f'{path}/Average_{pdb_chain}_repaired_BM.fxout'
            print(f'Processing file {fx_output_file}', file=sys.stderr)
            try:
                energy = pd.read_csv(fx_output_file, sep='\t', index_col=False, skiprows=8)

                with open(f'{path}/individual_list_0_PSSM.txt', 'r') as mut_file:
                    muts = mut_file.readlines()
            except FileNotFoundError as err:
                print(('WARNING: FoldX output or mutation list not found in FoldX directory '
                       '(contains .fxout files). Skipping'), file=sys.stderr)
                print(err, file=sys.stderr)
                continue

            if not len(energy) == len(muts):
                print(("WARNING: Number of mutations != number of predictions, suggesting "
                       "FoldX didn't complete correctly. Skipping"), file=sys.stderr)
                continue

            energy.columns = map(str.lower, energy.columns)
            energy = energy.drop('pdb', 1)
            energy['uniprot_id'] = uniprot_id
            energy['pdb_id'] = pdb_id
            energy['chain'] = chain

            energy['mut'] = muts
            energy['wt'] = energy.mut.str.get(0)
            energy['pos'] = energy.mut.str.slice(2, -3).astype(int)
            energy['mut'] = energy.mut.str.slice(-3, -2)

            energy = energy[['uniprot_id', 'pdb_id', 'chain', 'pos', 'wt', 'mut',
                             'sd', 'total energy', 'backbone hbond', 'sidechain hbond',
                             'van der waals', 'electrostatics', 'solvation polar',
                             'solvation hydrophobic', 'van der waals clashes',
                             'entropy sidechain', 'entropy mainchain', 'sloop_entropy',
                             'mloop_entropy', 'cis_bond', 'torsional clash',
                             'backbone clash', 'helix dipole', 'water bridge',
                             'disulfide', 'electrostatic kon', 'partial covalent bonds',
                             'energy ionisation', 'entropy complex']]

            energy.to_csv(sys.stdout, sep='\t', header=write_header, na_rep='NA',
                          index=False)
            write_header = False





def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dirs', metavar='D', nargs='+',
                        help="Directory(s) to crawl for FoldX output")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
