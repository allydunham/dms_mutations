#!/usr/bin/env python3
"""
Extract secondary structure from a pdb file via DSSP, and format as an SS8
file in the style of Porter5 (https://github.com/mircare/Porter5/).
PDB files are processed directly. Each PDB file should be followed by the chain to
process
DM files are used to find an associated PDB file, assuming PDBs are
in dm_dir/PDB_ID/PDB_ID.pdb
Output is to stdout for a single file, otherwise to
"""
import argparse
import sys
from os.path import basename, dirname

from Bio.PDB import PDBParser, DSSP

import deep_mut_tools as dmt

# Order of classification probability columns in Porter5 output
PROB_COLS = ['G', 'H', 'I', 'E', 'B', 'C', 'S', 'T']

def write_ss8(pdb_file, target_chain, pdb_id=None, pdb_offset=0, outfile=sys.stdout):
    """Run DSSP on input PDB, filter to specified chain and write to file"""
    if pdb_id is None:
        pdb_id = basename(pdb_file)

    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(pdb_id, pdb_file)
    dssp = DSSP(structure[0], pdb_file, dssp='mkdssp')

    print('#', 'AA', 'SS', *PROB_COLS, sep='\t', file=outfile)
    for key in dssp.keys():
        chain, (_, pos, _) = key
        if chain == target_chain:
            amino_acid, sec_struc = dssp[key][1:3]

            # Assume coil in Porter5 output means unstructured
            if sec_struc == '-':
                sec_struc = 'C'

            print(int(pos) + pdb_offset, amino_acid, sec_struc,
                  *make_ss_vector(sec_struc), sep='\t', file=outfile)

def make_ss_vector(ss_target, key=None):
    """Generate a one hot secondary structure vector vector"""
    if key is None:
        key = PROB_COLS

    return [1 if x == ss_target else 0 for x in key]

def unpack_pdb_id(pdb_id_str):
    """Unpack a coded PDB string from the dm file format"""
    spl = pdb_id_str.split(':')

    pdb_id, chain = spl[0:2]
    offset = int(spl[2]) if len(spl) > 2 else 0

    return pdb_id, chain, offset


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
                    pdb_id, chain, offset = unpack_pdb_id(pdb_id)
                    pdb_file = f'{dirname(dm_path)}/{pdb_id}/{pdb_id}.pdb'

                    print(f'-- {pdb_id} chain {chain}... ', file=sys.stderr)

                    if args.stdout:
                        print(f'\nSecondary Structures for chain {chain}, PDB {pdb_file}')
                        write_ss8(pdb_file, chain, pdb_id, offset)

                    else:
                        outfile = f'{dirname(dm_path)}/{pdb_id}/{pdb_id}.ss8'
                        with open(outfile, 'w') as outfile:
                            write_ss8(pdb_file, target_chain=chain, pdb_id=pdb_id,
                                      pdb_offset=offset, outfile=outfile)

            else:
                print(f'no PDB IDs listed\n')

    else: # PDB input
        # TODO allow multiple chains?
        if len(args.input) % 2 != 0:
            raise ValueError('PDB input files must be followed by a chain specification')

        pdb_input = iter(args.input)
        for path, chain in zip(pdb_input, pdb_input):
            if args.stdout:
                print(f'\nSecondary Structures for chain {chain}, PDB {path}')
                write_ss8(path, chain)

            else:
                outfile = f"{dirname(path)}/{basename(path).replace('.pdb', '.ss8')}"
                with open(outfile, 'w') as outfile:
                    write_ss8(path, chain, outfile=outfile)




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
