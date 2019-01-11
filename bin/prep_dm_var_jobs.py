#!/usr/bin/env python3
"""
Prepare files and LSF jobs to analyse deep mutagenesis data using a variety
of variant effect prediction tools.
"""
import sys
import os
import logging
from datetime import datetime
import argparse
import evcouplings.utils as ev
import deep_mut_tools as dm
import dmt
from nested_dicts import nested_merge
from smart_open import smart_open

ROOT_DIR = '/nfs/research1/beltrao/ally/'
EV_CONFIG_PATH = '/Users/ally/Projects/mutations/meta/base_evcouplings_config.txt'
ENV_HUMAN_DB = ROOT_DIR + 'databases/envision/human_predicted_combined_20170925.csv'
ENV_MOUSE_DB = ROOT_DIR + 'databases/envision/mouse_predicted_combined_20171004.csv'
ENV_YEAST_DB = ROOT_DIR + 'databases/envision/yeast_predicted_2017-03-12.csv'
UNIREF100 = ROOT_DIR + 'databases/uniprot/uniref100/uniref100_2019_1.fasta'
UNIREF90 = ROOT_DIR + 'databases/uniprot/uniref90/uniref90_2019_1.fasta'
UNIPROT = ROOT_DIR + 'databases/uniprot/uniprot/uniprot_2019_1.fasta'
ROTABASE_PATH = '/Users/ally/Projects/mutations/rotabase.txt'
LOG_ROOT = ROOT_DIR + 'logs/'

def main(args):
    """Main script"""
    # Prepare dirs and output file
    if args.name and args.name not in os.listdir(args.log):
        batch_name = args.name
    else:
        batch_name = f"dm_var_pred_{datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')}"

    log_dir = f"{args.log.rstrip('/')}/{batch_name}"
    #os.makedirs(log_dir)

    if args.out == '-':
        out_path = f'{batch_name}.sh'
    else:
        out_path = args.out

    # # Prepare for selected actions
    # if args.sift4g:
    #     pass

    # if args.envision:
    #     pass

    # if args.foldx:
    #     pass

    # if args.evcouplings:
    #     ev_config = nested_merge(ev.config.parse_config(args.ev_options),
    #                              ev.config.read_config_file(args.ev_default))

    # if args.polyphen2:
    #     pass

    with smart_open(out_path, 'w') as script_file:
        print("### Deep Mutagenesis Variant Prediction ###",
              f"# Batch name: {batch_name}",
              f"# Writing logs to: {log_dir}",
              "# Preparing jobs for:",
              f"#\tSIFT4G: {args.sift4g}",
              f"#\tEnvision: {args.envision}",
              f"#\tFoldX: {args.foldx}",
              f"#\tEVCouplings: {args.evcouplings}",
              f"#\tPolyphen2: {args.polyphen2}",
              sep='\n', end='\n\n', file=script_file)

        for dm_path in args.dm:
            print(f"### Jobs for {dm_path} ###", file=script_file)
            try:
                dm_dir = '/'.join(dm_path.split('/')[:-1])
                deep = dm.read_deep_mut(dm_path)
                tasker = dmt.DMTaskSelecter(deep)
                gene_name = deep.meta_data['gene_name']

                if args.sift4g:
                    print('# SIFT4G', file=script_file)
                    tasker.sift4g(path=dm_dir, dm_file=dm_path, overwrite=False)

                    command = ' '.join(["sift4g",
                                        f"-q {dm_dir}/{gene_name}.fa",
                                        f"-d {args.sift_db}",
                                        f"--subst {dm_dir}/{gene_name}.subst",
                                        f"--out {dm_dir}"])

                    job = ' '.join([f"bsub -o {log_dir}/sift4g.%J",
                                    f"-e {log_dir}/sift4g.%J.err",
                                    f'-M {args.sift_ram} -R "rusage[mem={args.sift_ram}]]"',
                                    f"'{command}'"])
                    print(job, file=script_file, end='\n\n')

                if args.envision:
                    print('# Envision', file=script_file)

                if args.foldx:
                    print('# FoldX', file=script_file)

                if args.evcouplings:
                    print('# EVCouplings', file=script_file)

                if args.polyphen2:
                    print('# Polyphen2', file=script_file)

            except Exception as err:
                raise err



def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dm', metavar='D', nargs='+', help="Input dm files")
    parser.add_argument('--out', '-o',
                        help='Output path ("-" to autogenerate or stdout by default)')

    jobs = parser.add_argument_group('Mutational effect jobs to Prepare')

    jobs.add_argument('--sift4g', '-s', action='store_true', help='Prepare Sift4g jobs')
    jobs.add_argument('--envision', '-e', action='store_true', help='Prepare Envision jobs')
    jobs.add_argument('--foldx', '-f', action='store_true', help='Prepare FoldX jobs')
    jobs.add_argument('--evcouplings', '-v', action='store_true', help='Prepare EVCouplings jobs')
    jobs.add_argument('--polyphen2', '-p', action='store_true', help='Prepare Polyphen2 jobs')

    sift = parser.add_argument_group('SIFT4G Options')
    sift.add_argument('--sift_db', default=UNIREF90, help='Base EVCouplings config file')
    sift.add_argument('--sift_ram', default=4000, type=int, help='RAM to allocate to sift4g')

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

    lsf = parser.add_argument_group('Technical LSF Options')
    lsf.add_argument('--log', '-l', help='Path to root logging directory', default=LOG_ROOT)
    lsf.add_argument('--name', '-n', help='Batch name')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
