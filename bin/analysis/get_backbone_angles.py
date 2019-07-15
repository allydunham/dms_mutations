#!/usr/bin/env python3
"""
Extract backbone angles from dm or PDB files
"""
import argparse
import sys
from math import pi
from os.path import basename, dirname
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils import seq1
import deep_mut_tools as dmt

RAD_FACTOR = 180/pi

def write_backbone_angles(pdb_file, chain='A', pdb_id=None, offset=0, outfile=sys.stdout):
    """
    Write Psi/Phi angles from a pdb file
    """
    if pdb_id is None:
        pdb_id = basename(pdb_file).replace('.pdb', '')

    pdb_parser = PDBParser()
    polypeptide_builder = PPBuilder()

    structure = pdb_parser.get_structure(pdb_id, pdb_file)

    polypeptides = polypeptide_builder.build_peptides(structure[0][chain])

    print('pdb_id', 'chain', 'aa', 'position', 'phi', 'psi', sep='\t', file=outfile)
    for peptide in polypeptides:
        angles = peptide.get_phi_psi_list()
        for residue, (phi, psi) in zip(peptide, angles):
            print(pdb_id, chain, seq1(residue.get_resname()),
                  residue.get_id()[1] + offset,
                  'NA' if phi is None else phi * RAD_FACTOR,
                  'NA' if psi is None else psi * RAD_FACTOR,
                  sep='\t', file=outfile)

def main(args):
    """Main script"""
    if args.dm:
        for dm_path in args.input:
            print(f'Processing {dm_path}...', file=sys.stderr)

            dm_header = dmt.read_deep_mut_header(dm_path)
            if 'pdb_id' in dm_header.keys() and dm_header['pdb_id']:
                pdbs_ids = dm_header['pdb_id']
                if isinstance(pdbs_ids, str):
                    pdbs_ids = [pdbs_ids]

                for pdb_id in pdbs_ids:
                    pdb_id, chain, offset = dmt.unpack_pdb_id(pdb_id)
                    pdb_file = f'{dirname(dm_path)}/{pdb_id}/{pdb_id}.pdb'

                    print(f'-- {pdb_id} chain {chain}... ', file=sys.stderr)

                    if args.stdout:
                        print(f'\nBackbone angles for chain {chain}, PDB {pdb_file}')
                        write_backbone_angles(pdb_file, chain, pdb_id, offset)

                    else:
                        outfile = f'{dirname(dm_path)}/{pdb_id}/{pdb_id}.bb_angles'
                        with open(outfile, 'w') as outfile:
                            write_backbone_angles(pdb_file, chain=chain, pdb_id=pdb_id,
                                                  offset=offset, outfile=outfile)

            else:
                print(f'no PDB IDs listed\n')

    else: # PDB input
        # TODO allow multiple chains?
        if len(args.input) % 2 != 0:
            raise ValueError('PDB input files must be followed by a chain specification')

        pdb_input = iter(args.input)
        for path, chain in zip(pdb_input, pdb_input):
            if args.stdout:
                print(f'\nBackbone angles for chain {chain}, PDB {path}')
                write_backbone_angles(path, chain)

            else:
                outfile = f"{dirname(path)}/{basename(path).replace('.pdb', '.bb_angles')}"
                with open(outfile, 'w') as outfile:
                    write_backbone_angles(path, chain, outfile=outfile)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', metavar='I', nargs='+',
                        help="Input PDB file(s)/dm file(s)")

    parser.add_argument('--dm', '-d', action='store_true', help="Treat input files as dm files")

    parser.add_argument('--stdout', '-s', action='store_true', help="Write to stdout")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
