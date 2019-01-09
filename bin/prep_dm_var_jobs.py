#!/usr/bin/env python3
"""
Prepare files and LSF jobs to analyse deep mutagenesis data with sift4g, envision, foldx, evcouplings and polyphen2
"""
import argparse
import fileinput

def main(args):
    """Main script"""

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dir', metavar='D', help="List of directories to process")

    parser.add_argument('--sift4g', '-s', action='store_true', help='Prepare Sift4g jobs')
    parser.add_argument('--envision', '-e', action='store_true', help='Prepare Envision jobs')
    parser.add_argument('--foldx', '-f', action='store_true', help='Prepare FoldX jobs')
    parser.add_argument('--evcouplings', '-v', action='store_true', help='Prepare EVCouplings jobs')
    parser.add_argument('--polyphen2', '-p', action='store_true', help='Prepare Polyphen2 jobs')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(vars(ARGS))
