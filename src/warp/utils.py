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


def apply_secular_correction(star_name, jd, rv, jd_ref=None):
    results = query_gaia(star_name)
    if results is None:
        print("[WARN] Gaia data not found or invalid, querying simbad.")
        results = query_simbad(star_name)
    if results is None:
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


def get_latest_pipeline(ins_name, pipelines):
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
        print(
            f"[WARN] No accepted pipelines found for instrument {ins_name}, keeping {pipelines}.")
        return pipelines
    accepted_ = [c for c in pipelines if accepted_pipelines[ins_name] in c]
    if len(accepted_) != 0:

        print(
            f"[INFO] Accepted pipeline for {ins_name}: {accepted_}")
        return accepted_

    else:
        print(
            f"[WARN] Accepted pipeline {accepted_pipelines[ins_name]} not found for instrument {ins_name}, keeping {pipelines}.")
        return pipelines
