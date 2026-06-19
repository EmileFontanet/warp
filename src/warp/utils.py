import numpy as np
from .simbad import query_simbad
from .gaia import query_gaia
import astropy.units as u


def secular_acceleration(pmra, pmdec, parallax):
    pm = np.sqrt(pmra**2 + pmdec**2)

    # distance in pc
    parallax_arcsec = parallax.to(u.arcsec)
    if parallax_arcsec.value <= 0:
        raise ValueError("Invalid parallax: must be > 0")

    d_pc = (1/parallax_arcsec.value) * u.pc

    # correct vt formula (numeric, then attach units)
    vt = 4.74047 * (pm/parallax).value * u.km/u.s
    # secular acceleration in m/s/yr
    dvdt = (vt**2 / d_pc).to(u.m/u.s/u.yr)

    return dvdt


def get_astrometric_data(star_name, astrometric_table=None, verbose=True):
    if astrometric_table is not None:
        if (verbose):
            print(
                'Using provided astrometric table to retrieve proper motion and parallax.')
        try:
            import pandas as pd
            if isinstance(astrometric_table, pd.DataFrame):
                row = None
                if star_name in astrometric_table.index:
                    row = astrometric_table.loc[star_name]
                else:
                    for col in ['MAIN_ID', 'ID', 'name', 'star_name', 'user_specified_id']:
                        if col in astrometric_table.columns:
                            clean_col = astrometric_table[col].astype(
                                str).str.replace(' ', '').str.lower()
                            clean_name = str(star_name).replace(
                                ' ', '').lower()
                            mask = clean_col == clean_name
                            if mask.sum() > 0:
                                row = astrometric_table[mask].iloc[0]
                                break
                if row is not None:
                    pmra = row.get('pmra', row.get('PMRA', None))
                    pmdec = row.get('pmdec', row.get('PMDEC', None))
                    plx = row.get('plx_value', row.get(
                        'parallax', row.get('PLX_VALUE', None)))

                    if pmra is not None and pmdec is not None and plx is not None:
                        if verbose:
                            print(
                                f"Found astrometric data for {star_name} in provided table.")
                        return {
                            "pmra": float(pmra) * u.mas/u.yr,
                            "pmdec": float(pmdec) * u.mas/u.yr,
                            "parallax": float(plx) * u.mas,
                            "source": 'simbad',
                        }
            else:
                # Try astropy table
                row = None
                for col in ['MAIN_ID', 'ID', 'name', 'star_name', 'id']:
                    if col in astrometric_table.colnames:
                        col_vals = astrometric_table[col]
                        if len(col_vals) > 0 and hasattr(col_vals[0], 'decode'):
                            col_vals = [v.decode('utf-8') for v in col_vals]
                        else:
                            col_vals = [str(v) for v in col_vals]

                        clean_col = [v.replace(' ', '').lower()
                                     for v in col_vals]
                        clean_name = str(star_name).replace(' ', '').lower()

                        if clean_name in clean_col:
                            idx = clean_col.index(clean_name)
                            row = astrometric_table[idx]
                            break
                if row is not None:
                    pmra_key = 'pmra' if 'pmra' in row.colnames else 'PMRA'
                    pmdec_key = 'pmdec' if 'pmdec' in row.colnames else 'PMDEC'
                    plx_key = 'plx_value' if 'plx_value' in row.colnames else 'parallax' if 'parallax' in row.colnames else 'PLX_VALUE'

                    if pmra_key in row.colnames and pmdec_key in row.colnames and plx_key in row.colnames:
                        pmra = row[pmra_key]
                        pmdec = row[pmdec_key]
                        plx = row[plx_key]
                        if not hasattr(pmra, 'unit') or pmra.unit is None:
                            pmra = pmra * u.mas/u.yr
                        if not hasattr(pmdec, 'unit') or pmdec.unit is None:
                            pmdec = pmdec * u.mas/u.yr
                        if not hasattr(plx, 'unit') or plx.unit is None:
                            plx = plx * u.mas
                        return {
                            "pmra": pmra,
                            "pmdec": pmdec,
                            "parallax": plx,
                            "source": 'simbad',
                        }
        except Exception as e:
            if verbose:
                print(
                    f"[WARN] Could not extract data from astrometric_table for {star_name}: {e}")

    try:
        results = query_gaia(star_name, verbose=verbose)
    except Exception as e:
        if verbose:
            print(f"[WARN] Gaia query failed: {e}, querying simbad.")
        results = None
    if results is None:
        if verbose:
            print("[WARN] Gaia data not found or invalid, querying simbad.")
        try:
            results = query_simbad(star_name, verbose=verbose)
        except Exception as e:
            if verbose:
                print(f"[WARN] Simbad query failed: {e}.")
            results = None
    return results


def apply_secular_correction(star_name, jd, rv, jd_ref=None, verbose=True, astrometric_table=None):
    if star_name is None:
        if verbose:
            print(
                "[WARN] No star name provided, skipping secular acceleration correction.")
        return rv

    results = get_astrometric_data(
        star_name, astrometric_table=astrometric_table, verbose=verbose)

    if results is None:
        if verbose:
            print("Could not retrieve astrometric data from Gaia or Simbad.")
        return rv
    else:
        if verbose:
            print(
                f"Retrieved astrometric data from {results['source']}: pmra={results['pmra']}, pmdec={results['pmdec']}, parallax={results['parallax']}")
    # simbad_data = query_simbad(star_name)
    if jd_ref is None and results['source'] == 'gaia':
        jd_ref = 57389.0
    elif jd_ref is None and results['source'] == 'simbad':
        jd_ref = 55500.0
    pmra = results['pmra']
    pmdec = results['pmdec']
    parallax = results['parallax']

    dvdt = secular_acceleration(pmra, pmdec, parallax)  # m/s/yr

    years = (jd - jd_ref) * (u.day.to(u.yr))

    correction = dvdt.value * years  # m/s

    return rv - correction


def get_latest_pipeline(ins_name, pipelines, verbose=True, skip_ndrs=True):
    from .config import accepted_pipelines
    from .config import coralie_ndrs_pipelines
    """
    Get the latest pipeline version for a given instrument.

    Args:
        ins_name (str): Instrument name.
        pipelines (list of str): List of pipeline versions. 
    Returns:
        str: Latest pipeline version for the instrument.
    """
    if ins_name not in accepted_pipelines:
        if verbose:
            print(
                f"[WARN] No accepted pipelines found for instrument {ins_name}, keeping {pipelines}.")
        return pipelines
    extra_drs = [] if skip_ndrs else coralie_ndrs_pipelines
    accepted_ = [
        c for c in pipelines if c in accepted_pipelines[ins_name] or c in extra_drs]

    # accepted_ = [c for c in pipelines if accepted_pipelines[ins_name] in c]
    if len(accepted_) != 0:
        if verbose:
            print(
                f"[INFO] Accepted pipeline for {ins_name}: {accepted_}")
        return accepted_

    else:
        if verbose:
            print(
                f"[WARN] Accepted pipeline {accepted_pipelines[ins_name]} not found for instrument {ins_name}, keeping {pipelines}.")
        return pipelines


def index_matching(list1, list2):
    '''
    Look for entries in list 2 that contain elements of list 1.

    :param list1: The list containing the substrings to search for
    :param list2: The list to search within
    '''
    if not isinstance(list1, list):
        list1 = list(list1)
    if not isinstance(list2, list):
        list2 = list(list2)
    return np.array([list2.index(c) for c in list2 if any(
        [f in c for f in list1])])


def doppler_shift(wave: np.ndarray, rv: float) -> np.ndarray:
    """
    Performs the doppler shift on the wavelength values.

    Args:
        wave (np.ndarray): the original wavelength values
        rv (float): the radial velocity of the object in km/s

    Returns:
        wave_shifted (np.ndarray): the doppler shifted wavelength values
    """
    from astropy.constants import c

    wave_shifted = wave + wave * rv / (c / 1e3).value
    return wave_shifted
