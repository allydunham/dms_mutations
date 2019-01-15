#!/usr/bin/env python3
"""
Prepare files and LSF jobs to analyse deep mutagenesis data using a variety
of variant effect prediction tools.
Generated script includes relative paths, so must be run from its directory
"""
#import sys
import os
import sys
import logging
import argparse
import configparser
from ftplib import FTP
from datetime import datetime

import evcouplings.utils as ev

import deep_mut_tools as dm
import dmt
from nested_dicts import nested_merge
from smart_open import smart_open

CONFIG = configparser.ConfigParser()
CONFIG.read(sys.path[0] + '/prep_dm_jobs_config.ini')

def main(args):
    """Main script"""
    # Prepare dirs and output file
    batch_name = (args.name if args.name and args.name not in os.listdir(args.log) else
                  f"dm_var_pred_{datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')}")

    log_dir = f"{args.log.rstrip('/')}/{batch_name}"

    if not CONFIG['misc'].getboolean('dummy_run'):
        os.makedirs(log_dir)

    # Initiate Log
    log_file = args.log_file or f'{log_dir}/log.txt'

    numeric_log_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError(f'Invalid log level: {args.log_level}')

    logging.basicConfig(filename=log_file, level=numeric_log_level,
                        format='%(asctime)s %(levelname)s: %(message)s')

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
        logging.debug('Envision databases available: %s', env_dbs)

    if args.foldx:
        ftp = FTP('ftp.ebi.ac.uk')
        ftp.login('anonymous')

    if args.evcouplings:
        ev_config = nested_merge(ev.config.read_config_file(args.ev_config),
                                 ev.config.parse_config(args.ev_options))
        logging.info('EVCouplings config file loaded')

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
                dm_id = dm_dir.split('/')[-1]

                logging.info('Loading DeepMut file: %s', dm_path)
                deep = dm.read_deep_mut(dm_path)

                logging.info('Initialising DMTaskSelecter for %s', dm_id)
                tasker = dmt.DMTaskSelecter(deep)

                logging.info('Writing Reference Fasta')
                tasker.ref_fasta(path=f"{dm_dir}/{deep.meta_data['gene_name']}.fa")



                if args.sift4g:
                    print('## SIFT4G', file=script_file)

                    logging.info('Preparing SIFT4G Files')
                    tasker.sift4g(path=dm_dir, dm_file=dm_path, overwrite=False)

                    logging.info('Writing SIFT4G Job')
                    print(sift_job(deep_mut=deep, out_dir=dm_dir, log_dir=log_dir,
                                   sift_db=args.sift_db, ram=args.sift_ram, dm_id=dm_id,
                                   batch_id=batch_name),
                          file=script_file, end='\n\n')

                if args.envision and deep.meta_data['species'] in env_dbs:
                    print('## Envision', file=script_file)

                    logging.info('Writing Envision Job')
                    print(envision_job(deep_mut=deep, dm_path=dm_path, out_dir=dm_dir,
                                       env_dbs=env_dbs, log_dir=log_dir, ram=args.env_ram,
                                       dm_id=dm_id, batch_id=batch_name),
                          file=script_file, end='\n\n')
                elif args.envision:
                    logging.warning('No envision database for %s', deep.meta_data['species'])

                if args.foldx and deep.meta_data['pdb_id']:
                    print('## FoldX', file=script_file)

                    logging.info('Preparing files for FoldX')
                    tasker.foldx(path=dm_dir, dm_file=dm_path, rotabase=args.rotabase, ftp=ftp)

                    for pdb in deep.meta_data['pdb_id']:
                        pdb = pdb.split(':')
                        logging.info('Writing FoldX jobs for %s', pdb[0])
                        print(f'# {pdb[0]}', file=script_file)
                        print(foldx_job(pdb_id=pdb[0], out_dir=dm_dir, log_dir=log_dir,
                                        ram=args.foldx_ram, dm_id=dm_id,
                                        batch_id=batch_name),
                              file=script_file, end='\n\n')
                elif args.foldx:
                    logging.warning('No PDB IDs in %s', dm_path)

                if args.evcouplings:
                    print('## EVCouplings', file=script_file)

                    logging.info('Preparing files for EVCouplings')
                    tasker.evcouplings(path=dm_dir, overwrite=False, ev_default=ev_config,
                                       ev_options='')

                    logging.info('Writing EVCouplings Job')
                    print(evcouplings_job(config=f'{dm_dir}/ev_config.txt', log_dir=log_dir,
                                          ram=args.ev_ram, dm_id=dm_id, batch_id=batch_name),
                          file=script_file, end='\n\n')

                if args.polyphen2:
                    print('## Polyphen2', file=script_file)

                    logging.info('Preparing files for Polyphen2')
                    tasker.polyphen2(path=f'{dm_dir}/polyphen2_variants.tsv')

                    logging.info('Writing jobs for Polyphen2')
                    print(*polyphen2_job(dm_dir=dm_dir, log_dir=log_dir, ram=args.pph_ram,
                                         gene_name=deep.meta_data['gene_name'], dm_id=dm_id,
                                         batch_id=batch_name),
                          file=script_file, sep='\n', end='\n\n')

            except Exception as err:
                logging.exception('Raised an exception while processing %s', dm_path)
                if args.foldx:
                    ftp.quit()
                raise err

    if args.foldx:
        ftp.quit()

def bsub(command, log, ram=8000, group=CONFIG['lsf']['lsf_group'], name='', dep='', cwd=''):
    """LSF submission command string"""
    command = ['bsub',
               f'-g {group}',
               f'-M {ram} -R "rusage[mem={ram}]"',
               f'-o {log}.%J',
               f"-e {log}.%J.err",
               f'"{command}"']

    if name:
        command.insert(1, f'-J "{name}"')

    if dep:
        command.insert(3, f'-w "{dep}"')

    if cwd:
        command.insert(4, f'-cwd "{cwd}"')

    return ' '.join(command)

def sift_job(deep_mut, out_dir, log_dir, sift_db, ram, dm_id, batch_id):
    """Generate LSF job string for SIFT4G"""
    command = ' '.join(["sift4g",
                        f"-q {out_dir}/{deep_mut.meta_data['gene_name']}.fa",
                        f"-d {sift_db}",
                        f"--subst {out_dir}/",
                        f"--out {out_dir}"])

    job_name = f'{batch_id}_{dm_id}_sift'

    return bsub(command, f'{log_dir}/{dm_id}_sift4g', ram=ram, name=job_name)

def envision_job(deep_mut, out_dir, dm_path, env_dbs, log_dir, ram, dm_id, batch_id):
    """Generate LSF job string for Envision"""
    uniprot_id = deep_mut.meta_data['uniprot_id']
    gene_name = deep_mut.meta_data['gene_name']
    species = deep_mut.meta_data['species']

    try:
        env_db = env_dbs[species]
    except KeyError:
        raise ValueError(f'No Envision Database for {species}')

    grep_command = ' '.join([CONFIG['misc']['table_grep'],
                             uniprot_id,
                             env_db,
                             f'> {out_dir}/{uniprot_id}_envision_db.csv'])

    py_command = ' '.join(['python',
                           CONFIG['misc']['dmt_path'],
                           f'--env {out_dir}/{uniprot_id}_envision_db.csv',
                           f'--path {out_dir}/{gene_name}_envision_vars.csv',
                           'envision',
                           dm_path])

    job_name = f'{batch_id}_{dm_id}_envision'

    return bsub(f'{grep_command};{py_command}', f'{log_dir}/{dm_id}_envision',
                ram=ram, name=job_name)

def foldx_job(pdb_id, out_dir, log_dir, ram, dm_id, batch_id):
    """Generate LSF job string for FoldX"""
    pdb_dir = f'{out_dir}/{pdb_id}'
    repair = ' '.join(['foldx',
                       '--command=RepairPDB',
                       f'--pdb={pdb_id}.pdb',
                       '--clean-mode=3'])

    model = ' '.join(['foldx',
                      '--command=BuildModel',
                      f'--pdb={pdb_id}_Repair.pdb',
                      f'--mutant-file=individual_list_{pdb_id}.txt',
                      '--numberOfRuns=3',
                      '--clean-mode=3',
                      '--out-pdb=false'])

    job_name = f'{batch_id}_{dm_id}_foldx'

    return bsub(f'{repair};{model}', f'{log_dir}/{dm_id}_foldx', ram=ram,
                name=job_name, cwd=pdb_dir)

def evcouplings_job(config, log_dir, ram, dm_id, batch_id):
    """Generate LSF job string for EVCouplings"""
    return bsub(f'evcouplings_runcfg {config}', f'{log_dir}/{dm_id}_evcouplings',
                ram=ram, name=f'{batch_id}_{dm_id}_evcouplings')

def polyphen2_job(dm_dir, gene_name, log_dir, ram, dm_id, batch_id):
    """Generate LSF job string for Polyphen2"""
    run_pph = ' '.join(['run_pph.pl',
                        f'{dm_dir}/polyphen2_variants.tsv',
                        f'1>{dm_dir}/pph_{gene_name}.features',
                        f'2>{dm_dir}/pph_{gene_name}.log'])

    run_weka = ' '.join(['run_weka.pl',
                         f'{dm_dir}/pph_{gene_name}.features',
                         f'1>{dm_dir}/pph_{gene_name}.predictions'])

    pph_name = f'{batch_id}_{dm_id}_polyphen2_pph'
    weka_name = f'{batch_id}_{dm_id}_polyphen2_weka'

    pph_job = bsub(run_pph, ram=ram,
                   log=f'{log_dir}/{dm_id}_polyphen2_pph',
                   name=pph_name)

    weka_job = bsub(run_weka, ram=ram,
                    log=f'{log_dir}/{dm_id}_polyphen2_weka',
                    name=weka_name, dep=f'ended({pph_name})')

    return [pph_job, weka_job]

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dm', metavar='D', nargs='+', help="Input dm files")
    parser.add_argument('--out', '-o',
                        help='Output path ("-" to autogenerate or stdout by default)')
    parser.add_argument('--log_file', help='Path to log file')
    parser.add_argument('--log_level', default='INFO', help='Path to log file')

    jobs = parser.add_argument_group('Mutational effect jobs to Prepare')

    jobs.add_argument('--sift4g', '-s', action='store_true', help='Prepare Sift4g jobs')
    jobs.add_argument('--envision', '-e', action='store_true', help='Prepare Envision jobs')
    jobs.add_argument('--foldx', '-f', action='store_true', help='Prepare FoldX jobs')
    jobs.add_argument('--evcouplings', '-v', action='store_true', help='Prepare EVCouplings jobs')
    jobs.add_argument('--polyphen2', '-p', action='store_true', help='Prepare Polyphen2 jobs')

    sift = parser.add_argument_group('SIFT4G Options')
    sift.add_argument('--sift_db', default=CONFIG['db']['uniref90'],
                      help='SIFT4G reference database')
    sift.add_argument('--sift_ram', default=4000, type=int, help='SIFT4G job RAM')

    evcoup = parser.add_argument_group('EVCouplings Options')
    evcoup.add_argument('--ev_config', default=CONFIG['evcouplings']['ev_config_path'],
                        help='Base EVCouplings config file')
    evcoup.add_argument('--ev_options', default='{}',
                        help='Extra EVCouplings options (overwrites ev_config parameters)')
    evcoup.add_argument('--ev_ram', default=8000, type=int, help='EVCouplings Job RAM')

    envision = parser.add_argument_group('Envision Options')
    envision.add_argument('--env_human', default=CONFIG['envision']['env_human_db'],
                          help='Envision human database path')
    envision.add_argument('--env_mouse', default=CONFIG['envision']['env_mouse_db'],
                          help='Envision mouse database path')
    envision.add_argument('--env_yeast', default=CONFIG['envision']['env_yeast_db'],
                          help='Envision yeast database path')
    envision.add_argument('--env_ram', default=1000, type=int, help='Envision job RAM')

    foldx = parser.add_argument_group('FoldX Options')
    foldx.add_argument('--rotabase', '-r', help='Path to FoldX rotabase.txt file',
                       default=CONFIG['foldx']['rotabase_path'])
    foldx.add_argument('--foldx_ram', default=8000, type=int, help='FoldX job RAM')

    pph = parser.add_argument_group('Polyphen2 Options')
    pph.add_argument('--pph_ram', default=1000, type=int, help='Polyphen2 job RAM')

    lsf = parser.add_argument_group('Technical LSF Options')
    lsf.add_argument('--log', '-l', help='Path to root logging directory',
                     default=CONFIG['lsf']['log_root'])
    lsf.add_argument('--group', help='LSF Group',
                     default=CONFIG['lsf']['lsf_group'])
    lsf.add_argument('--name', '-n', help='Batch name')

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
