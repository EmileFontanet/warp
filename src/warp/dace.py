import os
import argparse
import numpy as np
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from dace_query.spectroscopy import Spectroscopy
from astroquery.simbad import Simbad
import pandas as pd
from .simbad import query_simbad
from .utils import apply_secular_correction


def download_points(star, instrument=None, do_secular_corr=True,
                    skip_ndrs=True):
    """
    Download old and new DRS data for the specified star.
    Args:
        star_name (str): Name of the star to download data for.
    Returns:

    """
    if isinstance(instrument, str):
        instrument = [instrument]
    excluded_nights = ['2023-12-06',
                       '2023-12-02']  # List of nights to exclude, if any
    print(f"[INFO] Downloading data for star: {star.name}")
    dace_id = get_dace_id(star)
    filters = {
        'obj_id_daceid': {
            'equal': [dace_id]
        },
        'ins_name': {
            'contains': instrument
        }
    }
    if (skip_ndrs):
        filters['ins_name']['notContains'] = ['NDRS']

    results = Spectroscopy.query_database(
        filters=filters,
        output_format='pandas',
    )
    if results is None or len(results) == 0:
        raise ValueError(f"No data found for star {star.name} in DACE.")
    if (do_secular_corr):
        print("[INFO] Applying secular acceleration correction...")
        results['spectro_ccf_rv'] = apply_secular_correction(
            star.name,
            results['obj_date_bjd'],
            results['spectro_ccf_rv']
        )
    # results = results[~results.date_night.isin(excluded_nights)].copy()
    for ins in results.ins_name.unique():
        n_points = len(results[results.ins_name == ins])
        bad_qc_points = len(
            results[(results.ins_name == ins) & (~results.spectro_drs_qc)])
        print(
            f"[INFO] Retrieved {n_points} points, including {bad_qc_points} for which the QC failed, for instrument {ins}.")

    return results


def get_dace_id(star):
    result_table = Simbad.query_objectids(star.name)
    names = [star.name]
    try:
        hd_name = [id for id in result_table['id'] if id.startswith('HD')][0]
        star.hd = hd_name
        names.append(hd_name)
    except IndexError:
        star.hd = None
        print('[WARN] No HD name found for this star.')
    try:
        hip_name = [id for id in result_table['id'] if id.startswith('HIP')][0]
        star.hip = hip_name
        names.append(hip_name)
    except IndexError:
        star.hip = None
        print('[WARN] No HIP name found for this star.')
    all_names = names + list(result_table['id'].data)
    for name in all_names:
        name = name.replace(" ", "")
        print(
            f'[INFO] Trying to find DACE ID for star {star.name} with name {name}...')
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
            print(
                f'Could not find dace ID for star {star.name} with HD name {name}. Trying other name...')
    raise ValueError(
        f'Could not find DACE ID for star {star.name} with any known name.')


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
        print(
            f"[INFO] Performing MAD clipping for instrument: {ins} (threshold={threshold})")

        ins_data = data[data.ins_name == ins].copy()
        results[ins] = []

        for i in range(n_iter):
            median = np.median(ins_data.spectro_ccf_rv)
            mad = np.median(np.abs(ins_data.spectro_ccf_rv - median))
            if mad == 0:
                print(
                    f"[WARN] MAD = 0 for {ins} at iteration {i+1}, skipping.")
                break

            modified_z = 0.6745 * (ins_data.spectro_ccf_rv - median) / mad
            mask = np.abs(modified_z) < threshold

            filtered_data = ins_data[mask]
            removed_data = ins_data[~mask]

            print(
                f"[INFO] Iteration {i+1}: clipped {len(removed_data)} outliers.")

            results[ins].append((filtered_data.copy(), removed_data.copy()))

            # If no outliers were removed, stop early
            if len(removed_data) == 0:
                break

            ins_data = filtered_data.copy()

    return results
