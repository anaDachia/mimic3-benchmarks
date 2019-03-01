from __future__ import absolute_import
from __future__ import print_function

import argparse

import os
import sys

from mimic3benchmark.subject import read_stays, get_events_for_stay, add_hours_elpased_to_events, read_seq_table
from mimic3benchmark.subject import convert_seq_tables_to_timeseries


parser = argparse.ArgumentParser(description='Extract seq episodes from per-subject data.')
parser.add_argument('subjects_root_path', type=str, help='Directory containing subject sub-directories.')
parser.add_argument('--seq_tables', '-s', type=str, nargs='+', help='Tables from which to read sequential tables.',
                    default=['LABEVENTS', 'MICROBIOLOGYEVENTS']) #, 'PRESCRIPTION'
parser.add_argument('--verbose', '-v', type=int, help='Level of verbosity in output.', default=1)
args, _ = parser.parse_known_args()


for subject_dir in os.listdir(args.subjects_root_path):
    dn = os.path.join(args.subjects_root_path, subject_dir)
    try:
        subject_id = int(subject_dir)
        if not os.path.isdir(dn):
            raise Exception
    except:
        continue
    sys.stdout.write('Subject {}: '.format(subject_id))
    sys.stdout.flush()
    try:
        sys.stdout.write('reading...')
        sys.stdout.flush()
        stays = read_stays(os.path.join(args.subjects_root_path, subject_dir))
    except:
        sys.stdout.write('error reading from disk!\n')
        continue

    sys.stdout.write('extracting separate episodes...')
    sys.stdout.flush()
    for i in range(stays.shape[0]):
        stay_id = stays.ICUSTAY_ID.iloc[i]
        sys.stdout.write(' {}'.format(stay_id))
        sys.stdout.flush()
        intime = stays.INTIME.iloc[i]
        outtime = stays.OUTTIME.iloc[i]

        for seq_table_name in args.seq_tables:
            if not os.path.exists(os.path.join(args.subjects_root_path, subject_dir, seq_table_name.lower() + '.csv')):
                print ("no event of this type for this subject: " + str(seq_table_name))
                continue
            seq_table = read_seq_table(os.path.join(args.subjects_root_path, subject_dir), seq_table_name)
            timeseries = convert_seq_tables_to_timeseries(seq_table)
            episode = get_events_for_stay(timeseries, stay_id, intime, outtime)
            if episode.shape[0] == 0:
                sys.stdout.write(' (no data!)')
                sys.stdout.flush()
                continue
            episode = add_hours_elpased_to_events(episode, intime)#.set_index('HOURS').sort_index(axis=0)
            episode = episode.pivot(index= 'ITEMID', columns='HOURS', values='FLAG')


            columns = list(episode.columns)
            columns_sorted = sorted(columns)
            episode = episode[columns_sorted]
            episode.to_csv(os.path.join(args.subjects_root_path, subject_dir,
                                        'episode{}'.format(i+1), seq_table_name.lower() +'_timeseries.csv'), index_label='Hours')
    sys.stdout.write(' DONE!\n')

