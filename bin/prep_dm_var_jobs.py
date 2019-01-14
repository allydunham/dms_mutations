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

ROOT_DIR = '/nfs/research1/beltrao/ally'
EV_CONFIG_PATH = '/Users/ally/Projects/mutations/meta/base_evcouplings_config.txt'

ENV_HUMAN_DB = ROOT_DIR + '/databases/envision/human_predicted_combined_20170925.csv'
ENV_MOUSE_DB = ROOT_DIR + '/databases/envision/mouse_predicted_combined_20171004.csv'
ENV_YEAST_DB = ROOT_DIR + '/databases/envision/yeast_predicted_2017-03-12.csv'

UNIREF100 = ROOT_DIR + '/databases/uniprot/uniref100/uniref100_2019_1.fasta'
UNIREF90 = ROOT_DIR + '/databases/uniprot/uniref90/uniref90_2019_1.fasta'
UNIPROT = ROOT_DIR + '/databases/uniprot/uniprot/uniprot_2019_1.fasta'
ROTABASE_PATH = '/Users/ally/Projects/mutations/rotabase.txt'
LOG_ROOT = ROOT_DIR + '/logs'

def main(args):
    """Main script"""
    # Prepare dirs and output file
    batch_name = (args.name if args.name and args.name not in os.listdir(args.log) else
                  f"dm_var_pred_{datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')}")

    log_dir = f"{args.log.rstrip('/')}/{batch_name}"
    #os.makedirs(log_dir)

    # Initiate Log
    log_file = args.log_file if args.log_file else f'{log_dir}/log.txt'
    logging.basicConfig(filename=log_file, level=logging.INFO)
    logging.info('Batch name %s', batch_name)
    logging.info('Writing job logs to %s', log_dir)

    out_path = f'{batch_name}.sh' if args.out == '-' else args.out
    logging.info('Writing job script to %s', out_path if args.out else 'STDOUT')

    # Prepare for selected actions
    # if args.sift4g:
    #     pass

    if args.envision:
        env_dbs = {'Homo sapiens': args.env_human,
                   'Saccharomyces cerevisiae': args.env_yeast,
                   'Mus musculus': args.env_mouse}

    # if args.foldx:
    #     pass

    if args.evcouplings:
        ev_config = nested_merge(ev.config.read_config_file(args.ev_config),
                                 ev.config.parse_config(args.ev_options))

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
                tasker.ref_fasta(path=f"{dm_dir}/{deep.meta_data['gene_name']}.fa")

                if args.sift4g:
                    print('## SIFT4G', file=script_file)
                    tasker.sift4g(path=dm_dir, dm_file=dm_path, overwrite=False)
                    print(sift_job(deep_mut=deep, out_dir=dm_dir, log_dir=log_dir,
                                   sift_db=args.sift_db, ram=args.sift_ram),
                          file=script_file, end='\n\n')

                if args.envision and deep.meta_data['species'] in env_dbs:
                    print('## Envision', file=script_file)
                    print(envision_job(deep_mut=deep, dm_path=dm_path, out_dir=dm_dir,
                                       env_dbs=env_dbs, log_dir=log_dir, ram=args.env_ram),
                          file=script_file, end='\n\n')
                elif args.envision:
                    logging.warning('No envision database for %s', deep.meta_data['species'])

                if args.foldx and deep.meta_data['pdb_id']:
                    print('## FoldX', file=script_file)
                    tasker.foldx(path=dm_dir, dm_file=dm_path, rotabase=args.rotabase)
                    for pdb in deep.meta_data['pdb_id']:
                        pdb = pdb.split(':')
                        print(f'# {pdb[0]}', file=script_file)
                        print(foldx_job(pdb_id=pdb[0], out_dir=dm_dir, log_dir=log_dir,
                                        ram=args.foldx_ram),
                              file=script_file, end='\n\n')
                elif args.foldx:
                    logging.warning('No PDB IDs in %s', dm_path)

                if args.evcouplings:
                    print('## EVCouplings', file=script_file)
                    tasker.evcouplings(path=dm_dir, overwrite=False, ev_default=ev_config,
                                       ev_options='')
                    print(evcouplings_job(config=f'{dm_dir}/ev_config.txt', log_dir=log_dir,
                                          ram=args.ev_ram),
                          file=script_file, end='\n\n')

                if args.polyphen2:
                    print('## Polyphen2', file=script_file)

            except Exception as err:
                raise err

def sift_job(deep_mut, out_dir, log_dir, sift_db, ram):
    """Generate LSF job string for SIFT4G"""
    command = ' '.join(["sift4g",
                        f"-q {out_dir}/{deep_mut.meta_data['gene_name']}.fa",
                        f"-d {sift_db}",
                        f"--subst {out_dir}/{deep_mut.meta_data['gene_name']}.subst",
                        f"--out {out_dir}"])

    return' '.join([f"bsub -o {log_dir}/sift4g.%J",
                    f"-e {log_dir}/sift4g.%J.err",
                    f'-M {ram} -R "rusage[mem={ram}]]"',
                    f"'{command}'"])

def envision_job(deep_mut, out_dir, dm_path, env_dbs, log_dir, ram):
    """Generate LSF job string for Envision"""
    uniprot_id = deep_mut.meta_data['uniprot_id']
    gene_name = deep_mut.meta_data['gene_name']
    species = deep_mut.meta_data['species']

    try:
        env_db = env_dbs[species]
    except KeyError:
        raise ValueError(f'No Envision Database for {species}')



    grep_command = ' '.join(['cat',
                             f'<(head -n 1 {env_db})',
                             f'<(grep {uniprot_id} {env_db})',
                             f'> {out_dir}/{uniprot_id}_envision_db.csv'])

    py_command = ' '.join(['python',
                           f'{ROOT_DIR}/mutations/bin/dmt.py',
                           f'--env {out_dir}/{uniprot_id}_envision_db.csv',
                           f'--path {out_dir}/{gene_name}_envision_vars.csv',
                           'envision',
                           dm_path])

    return ' '.join([f"bsub -o {log_dir}/envision.%J",
                     f"-e {log_dir}/envision.%J.err",
                     f'-M {ram} -R "rusage[mem={ram}]]"',
                     f"'{grep_command};{py_command}'"])

def foldx_job(pdb_id, out_dir, log_dir, ram):
    """Generate LSF job string for FoldX"""
    pdb_dir = f'{out_dir}/{pdb_id}'
    repair = ' '.join(['foldx',
                       '--command=RepairPDB',
                       f'--pdb={pdb_dir}/{pdb_id}.pdb',
                       '--clean-mode=3'])

    model = ' '.join(['foldx',
                      '--command=BuildModel',
                      f'--pdb={pdb_dir}/{pdb_id}_Repair.pdb',
                      f'--mutant-file=individual_list_{pdb_id}.txt',
                      '--numberOfRuns=3',
                      '--clean-mode=3'])

    return ' '.join([f"bsub -o {log_dir}/envision.%J",
                     f"-e {log_dir}/envision.%J.err",
                     f'-M {ram} -R "rusage[mem={ram}]]"',
                     f"'{repair};{model}'"])

def evcouplings_job(config, log_dir, ram):
    """Generate LSF job string for EVCouplings"""
    command = f'evcouplings_runcfg {config}'

    return' '.join([f"bsub -o {log_dir}/sift4g.%J",
                    f"-e {log_dir}/sift4g.%J.err",
                    f'-M {ram} -R "rusage[mem={ram}]]"',
                    f"'{command}'"])

def polyphen2_job():
    """Generate LSF job string for Polyphen2"""


def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dm', metavar='D', nargs='+', help="Input dm files")
    parser.add_argument('--out', '-o',
                        help='Output path ("-" to autogenerate or stdout by default)')
    parser.add_argument('--log_file', help='Path to log file')

    jobs = parser.add_argument_group('Mutational effect jobs to Prepare')

    jobs.add_argument('--sift4g', '-s', action='store_true', help='Prepare Sift4g jobs')
    jobs.add_argument('--envision', '-e', action='store_true', help='Prepare Envision jobs')
    jobs.add_argument('--foldx', '-f', action='store_true', help='Prepare FoldX jobs')
    jobs.add_argument('--evcouplings', '-v', action='store_true', help='Prepare EVCouplings jobs')
    jobs.add_argument('--polyphen2', '-p', action='store_true', help='Prepare Polyphen2 jobs')

    sift = parser.add_argument_group('SIFT4G Options')
    sift.add_argument('--sift_db', default=UNIREF90, help='SIFT4G reference database')
    sift.add_argument('--sift_ram', default=4000, type=int, help='SIFT4G job RAM')

    evcoup = parser.add_argument_group('EVCouplings Options')
    evcoup.add_argument('--ev_config', default=EV_CONFIG_PATH,
                        help='Base EVCouplings config file')
    evcoup.add_argument('--ev_options', default='{}',
                        help='Extra EVCouplings options (overwrites ev_config parameters)')
    evcoup.add_argument('--ev_ram', default=8000, type=int, help='EVCouplings Job RAM')

    envision = parser.add_argument_group('Envision Options')
    envision.add_argument('--env_human', default=ENV_HUMAN_DB, help='Envision human database path')
    envision.add_argument('--env_mouse', default=ENV_MOUSE_DB, help='Envision mouse database path')
    envision.add_argument('--env_yeast', default=ENV_YEAST_DB, help='Envision yeast database path')
    envision.add_argument('--env_ram', default=1000, type=int, help='Envision job RAM')

    foldx = parser.add_argument_group('FoldX Options')
    foldx.add_argument('--rotabase', '-r', help='Path to FoldX rotabase.txt file',
                       default=ROTABASE_PATH)
    foldx.add_argument('--foldx_ram', default=8000, type=int, help='FoldX job RAM')

    lsf = parser.add_argument_group('Technical LSF Options')
    lsf.add_argument('--log', '-l', help='Path to root logging directory', default=LOG_ROOT)
    lsf.add_argument('--name', '-n', help='Batch name')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
