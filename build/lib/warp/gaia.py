
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import logging


def query_gaia(name, radius=2.0*u.arcsec, verbose=True):
    from astroquery.gaia import Gaia

    coord = SkyCoord(name) if isinstance(
        name, str) and ',' in name else SkyCoord.from_name(name)
    j = Gaia.cone_search(coord, radius=radius).get_results()
    if verbose is False:
        logging.getLogger('astroquery.utils.tap.core').setLevel(logging.ERROR)
        logging.getLogger('astroquery.utils.tap').setLevel(logging.ERROR)
        logging.getLogger('astroquery').setLevel(logging.ERROR)
    if len(j) > 0:
        # pick best (smallest angular separation / best RUWE)
        row = j[0]
        parallax = row['parallax'] * u.mas
        pmra = row['pmra'] * u.mas/u.yr
        pmdec = row['pmdec'] * u.mas/u.yr
        ruwe = row.get('ruwe', None)
        par_snr = (row['parallax'] / row['parallax_error']
                   ) if row['parallax_error'] else 0
        if (parallax.value > 0) and (par_snr >= 3) and (ruwe is None or ruwe < 1.4):
            return dict(source='gaia', parallax=parallax, pmra=pmra, pmdec=pmdec, ruwe=ruwe)
    return None


def query_gaia_photometry(gaia_id):
    from astroquery.gaia import Gaia
    job = Gaia.launch_job("""
    SELECT source_id, ra, dec
    FROM gaiadr3.gaia_source
    WHERE source_id IN (...)
    """)
    sources = job.get_results()


def get_gaia_dr3_radii(source_ids, chunk_size=5000):
    from astroquery.gaia import Gaia
    """
    Query Gaia DR3 GSP-Phot radii for a list of Gaia DR3 source_ids.

    Parameters
    ----------
    df : pandas.DataFrame
        Input table containing Gaia DR3 source IDs.
    source_id_col : str
        Name of the source_id column.
    chunk_size : int
        Number of stars queried at once.

    Returns
    -------
    pandas.DataFrame
        Original dataframe merged with Gaia DR3 radius and Teff columns.
    """

    source_ids = (
        source_ids
        .dropna()
        .astype("int64")
        .astype(str)
        .unique()
    )

    results = []

    for i in range(0, len(source_ids), chunk_size):
        chunk = source_ids[i:i + chunk_size]
        ids = ",".join(chunk)

        query = f"""
        SELECT
            source_id,
            radius_gspphot,
            radius_gspphot_lower,
            radius_gspphot_upper,
            teff_gspphot,
            teff_gspphot_lower,
            teff_gspphot_upper
        FROM gaiadr3.astrophysical_parameters
        WHERE source_id IN ({ids})
        """

        job = Gaia.launch_job_async(query)
        tab = job.get_results().to_pandas()
        results.append(tab)

    if len(results) == 0:
        return None

    gaia = pd.concat(results, ignore_index=True)
    gaia["source_id"] = gaia["source_id"].astype("int64")

    return gaia
