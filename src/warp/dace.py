import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from dace_query.spectroscopy import Spectroscopy
from astroquery.simbad import Simbad
import pandas as pd

def download_points(star_name):
    """
    Download old and new DRS data for the specified star.
    Args:
        star_name (str): Name of the star to download data for.
    Returns:
        
    """
    excluded_nights = ['2023-12-06', '2023-12-02']  # List of nights to exclude, if any
    print(f"[INFO] Downloading data for star: {star_name}")
    dace_id = get_dace_id(star_name)
    filter = {
        'obj_id_daceid': {
            'equal': [dace_id]
        },
        'ins_name' : {
            'contains' : ['CORALIE14', 'CORALIE24']
        }
    }
    results = pd.DataFrame.from_dict(Spectroscopy.query_database(filters=filter, limit=1000))
    results = results[~results.date_night.isin(excluded_nights)].copy()
    for ins in results.ins_name.unique():
        n_points = len(results[results.ins_name == ins])
        bad_qc_points = len(results[(results.ins_name == ins) & (~results.spectro_drs_qc)])
        print(f"[INFO] Retrieved {n_points} points, including {bad_qc_points} for which the QC failed, for instrument {ins}.")
    filtered_results = results[(results.spectro_drs_qc) & (results.spectro_ccf_rv_err > 0)].copy()
    return results, filtered_results

def get_dace_id(star_name):
    result_table = Simbad.query_objectids(star_name)
    try:
        hd_name = [id for id in result_table['id'] if id.startswith('HD')][0]
    except IndexError:
        hd_name = None
        print('[WARN] No HD name found for this star.')
    try:
        hip_name = [id for id in result_table['id'] if id.startswith('HIP')][0]
    except IndexError:
        hip_name = None
        print('[WARN] No HIP name found for this star.')
    for name in [star_name, hd_name, hip_name]:
        name = name.replace(" ", "")
        print(f'[INFO] Trying to find DACE ID for star {star_name} with name {name}...')
        try:
            filter = {
                'obj_id_catname': {
                    'equal': [name]
                }
            }
            results = Spectroscopy.query_database(filters=filter, limit=10)
            dace_id = results['obj_id_daceid'][results['obj_id_catname'] == name][0]
            return dace_id
        except Exception as e:
            print(f'Could not find dace ID for star {star_name} with HD name {name}. Trying other name...')
    raise ValueError(f'Could not find DACE ID for star {star_name} with any known name.')
    

def mad_clip(data, threshold=5, n_iter=3):
    """
    Perform Median Absolute Deviation (MAD) clipping iteratively for each instrument.
    Args:
        data (pd.DataFrame): DataFrame with at least ['ins_name', 'spectro_ccf_rv'] columns.
        threshold (float): MAD clipping threshold.
        n_iter (int): Number of iterations per instrument.
    Returns:
        dict: Dictionary with filtered data and removed outliers at each iteration.
              { instrument: [ (filtered_data, removed_data), ... ] }
    """
    results = {}
    data = data.copy()

    for ins in data.ins_name.unique():
        print(f"[INFO] Performing MAD clipping for instrument: {ins} (threshold={threshold})")

        ins_data = data[data.ins_name == ins].copy()
        results[ins] = []

        for i in range(n_iter):
            median = np.median(ins_data.spectro_ccf_rv)
            mad = np.median(np.abs(ins_data.spectro_ccf_rv - median))
            if mad == 0:
                print(f"[WARN] MAD = 0 for {ins} at iteration {i+1}, skipping.")
                break

            modified_z = 0.6745 * (ins_data.spectro_ccf_rv - median) / mad
            mask = np.abs(modified_z) < threshold

            filtered_data = ins_data[mask]
            removed_data = ins_data[~mask]

            print(f"[INFO] Iteration {i+1}: clipped {len(removed_data)} outliers.")

            results[ins].append((filtered_data.copy(), removed_data.copy()))

            # If no outliers were removed, stop early
            if len(removed_data) == 0:
                break

            ins_data = filtered_data.copy()

    return results