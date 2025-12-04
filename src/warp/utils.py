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


def apply_secular_correction(star_name, jd, rv, jd_ref=None, verbose=True):
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
    if results is None:
        if verbose:
            print("Could not retrieve astrometric data from Gaia or Simbad.")
        return rv
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


def get_latest_pipeline(ins_name, pipelines, verbose=True):
    from .config import accepted_pipelines
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
    accepted_ = [c for c in pipelines if accepted_pipelines[ins_name] in c]
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
