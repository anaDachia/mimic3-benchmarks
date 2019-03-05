from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import os
import pandas as pd

from mimic3benchmark.util import *


def read_stays(subject_path):
    stays = dataframe_from_csv(os.path.join(subject_path, 'stays.csv'), index_col=None)
    stays.INTIME = pd.to_datetime(stays.INTIME)
    stays.OUTTIME = pd.to_datetime(stays.OUTTIME)
    stays.DOB = pd.to_datetime(stays.DOB)
    stays.DOD = pd.to_datetime(stays.DOD)
    stays.DEATHTIME = pd.to_datetime(stays.DEATHTIME)
    stays.sort_values(by=['INTIME', 'OUTTIME'], inplace=True)
    return stays


def read_diagnoses(subject_path):
    return dataframe_from_csv(os.path.join(subject_path, 'diagnoses.csv'), index_col=None)


def read_events(subject_path, remove_null=True):
    events = dataframe_from_csv(os.path.join(subject_path, 'events.csv'), index_col=None)
    if remove_null:
        events = events.ix[events.VALUE.notnull()]
    events.CHARTTIME = pd.to_datetime(events.CHARTTIME)
    events.HADM_ID = events.HADM_ID.fillna(value=-1).astype(int)
    events.ICUSTAY_ID = events.ICUSTAY_ID.fillna(value=-1).astype(int)
    events.VALUEUOM = events.VALUEUOM.fillna('').astype(str)
    # events.sort_values(by=['CHARTTIME', 'ITEMID', 'ICUSTAY_ID'], inplace=True)
    return events

def read_seq_table(subject_path, table_name):
    events = dataframe_from_csv(os.path.join(subject_path, table_name + '.csv'), index_col=None)
    events.CHARTTIME = pd.to_datetime(events.CHARTTIME)
    events.HADM_ID = events.HADM_ID.fillna(value=-1).astype(int)
    events.ICUSTAY_ID = events.ICUSTAY_ID.fillna(value=-1).astype(int)
    return events

def read_static_table(subject_path, table_name):
    events = dataframe_from_csv(os.path.join(subject_path, table_name.lower() + '.csv'), index_col=None)
    events.HADM_ID = events.HADM_ID.fillna(value=-1).astype(int)
    events.SUBJECT_ID = events.SUBJECT_ID.fillna(value=-1).astype(int)
    return events


def get_events_for_stay(events, icustayid, intime=None, outtime=None):
    idx = (events.ICUSTAY_ID == icustayid)
    if intime is not None and outtime is not None:
        idx = idx | ((events.CHARTTIME >= intime) & (events.CHARTTIME <= outtime))
    events = events.ix[idx]
    del events['ICUSTAY_ID']
    return events


def add_hours_elpased_to_events(events, dt, remove_charttime=True):
    events['HOURS'] = (events.CHARTTIME - dt).apply(lambda s: s / np.timedelta64(1, 's')) / 60./60
    if remove_charttime:
        del events['CHARTTIME']
    return events


def convert_events_to_timeseries(events, variable_column='VARIABLE', variables=[]):
    metadata = events[['CHARTTIME', 'ICUSTAY_ID']].sort_values(by=['CHARTTIME', 'ICUSTAY_ID'])\
                    .drop_duplicates(keep='first').set_index('CHARTTIME')
    timeseries = events[['CHARTTIME', variable_column, 'VALUE']]\
                    .sort_values(by=['CHARTTIME', variable_column, 'VALUE'], axis=0)\
                    .drop_duplicates(subset=['CHARTTIME', variable_column], keep='last')
    timeseries = timeseries.pivot(index='CHARTTIME', columns=variable_column, values='VALUE').merge(metadata, left_index=True, right_index=True)\
                    .sort_index(axis=0).reset_index()
    for v in variables:
        if v not in timeseries:
            timeseries[v] = np.nan
    return timeseries

def convert_seq_tables_to_timeseries(seq_table_events):
    #print(list(seq_table_events))
    metadata = seq_table_events[['CHARTTIME', 'ICUSTAY_ID']].sort_values(by=['CHARTTIME', 'ICUSTAY_ID'])\
                    .drop_duplicates(keep='first')#.set_index('CHARTTIME')
    #print(list(metadata))
    timeseries = seq_table_events[['CHARTTIME', "ITEMID", 'FLAG']]\
            .sort_values(by=['CHARTTIME', 'ITEMID', 'FLAG'], axis=0)\
            .drop_duplicates(subset=['CHARTTIME', 'ITEMID', 'FLAG'], keep='last')
    # timeseries = timeseries.pivot(index='CHARTTIME', columns='ITEMID', values='FLAG').merge(metadata, left_index=True, right_index=True)\
    #                 .sort_index(axis=0).reset_index()
    #print(list(timeseries))
    timeseries = timeseries.merge(metadata, on = "CHARTTIME").reset_index()
    return timeseries




def get_first_valid_from_timeseries(timeseries, variable):
    if variable in timeseries:
        idx = timeseries[variable].notnull()
        if idx.any():
            loc = np.where(idx)[0][0]
            return timeseries[variable].iloc[loc]
    return np.nan
