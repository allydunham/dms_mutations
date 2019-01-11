#!/usr/bin/env python3
"""
Prepare files and LSF jobs to analyse deep mutagenesis data using a variety
of variant effect prediction tools.
"""
import argparse
import evcouplings.utils as ev
import deep_mut_tools as dm
import dmt
from nested_dicts import nested_merge

ROOT_DIR = '/nfs/research1/beltrao/ally/'
EV_CONFIG_PATH = '/Users/ally/Projects/mutations/meta/base_evcouplings_config.txt'
ENV_HUMAN_DB = ROOT_DIR + 'databases/envision/human_predicted_combined_20170925.csv'
ENV_MOUSE_DB = ROOT_DIR + 'databases/envision/mouse_predicted_combined_20171004.csv'
ENV_YEAST_DB = ROOT_DIR + 'databases/envision/yeast_predicted_2017-03-12.csv'
ROTABASE_PATH = '/Users/ally/Projects/mutations/rotabase.txt'

def main(args):
    """Main script"""
    # Prepare for selected actions
    if args.sift4g:
        pass

    if args.envision:
        pass

    if args.foldx:
        pass

    if args.evcouplings:
        ev_config = nested_merge(ev.config.parse_config(args.ev_options),
                                 ev.config.read_config_file(args.ev_default))

    if args.polyphen2:
        pass

    # Parse folders and process .dm files found




def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dm', metavar='D', nargs='+', help="Input dm files")

    jobs = parser.add_argument_group('Mutational effect jobs to Prepare')

    jobs.add_argument('--sift4g', '-s', action='store_true', help='Prepare Sift4g jobs')
    jobs.add_argument('--envision', '-e', action='store_true', help='Prepare Envision jobs')
    jobs.add_argument('--foldx', '-f', action='store_true', help='Prepare FoldX jobs')
    jobs.add_argument('--evcouplings', '-v', action='store_true', help='Prepare EVCouplings jobs')
    jobs.add_argument('--polyphen2', '-p', action='store_true', help='Prepare Polyphen2 jobs')

    evcoup = parser.add_argument_group('EVCouplings Options')
    evcoup.add_argument('--ev_config', default=EV_CONFIG_PATH, help='Base EVCouplings config file')
    evcoup.add_argument('--ev_options', help='Additional EVCouplings config options')

    envision = parser.add_argument_group('Envision Options')
    envision.add_argument('--env_human', default=ENV_HUMAN_DB, help='Envision human database path')
    envision.add_argument('--env_mouse', default=ENV_MOUSE_DB, help='Envision mouse database path')
    envision.add_argument('--env_yeast', default=ENV_YEAST_DB, help='Envision yeast database path')

    foldx = parser.add_argument_group('FoldX Options')
    foldx.add_argument('--rotabase', '-r', help='Path to FoldX rotabase.txt file',
                        default=ROTABASE_PATH)

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
