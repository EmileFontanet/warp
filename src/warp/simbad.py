# src/warp/simbad.py
from astropy.coordinates import SkyCoord
import astropy.units as u
import logging


def query_simbad(name, verbose=True):
    from astroquery.simbad import Simbad

    """Return coordinates and basic astrometric params from SIMBAD."""
    if verbose is False:
        logging.getLogger('astroquery.utils.tap.core').setLevel(logging.ERROR)
        logging.getLogger('astroquery.utils.tap').setLevel(logging.ERROR)
        logging.getLogger('astroquery').setLevel(logging.ERROR)
    Simbad.reset_votable_fields()
    Simbad.add_votable_fields('pmra', 'pmdec', 'ra',
                              'dec', 'plx_value', 'rvz_radvel')

    result = Simbad.query_object(name)
    if result is None:
        raise ValueError(f"SIMBAD could not find object: {name}")

    # ra(d) and dec(d) are both in degrees
    ra = result["ra"][0]
    dec = result["dec"][0]
    coord = SkyCoord(ra, dec, unit=(u.deg, u.deg))

    pmra = result["pmra"][0] * u.mas/u.yr
    pmdec = result["pmdec"][0] * u.mas/u.yr
    plx = result["plx_value"][0] * u.mas
    rv = result["rvz_radvel"][0] * u.km/u.s
    return {
        "coord": coord,
        "pmra": pmra,
        "pmdec": pmdec,
        "parallax": plx,
        "rv": rv,
        "source": 'simbad'
    }


def get_ids(star):
    from astroquery.simbad import Simbad

    result_table = Simbad.query_objectids(star)
    return [result_table['id'][i] for i in range(len(result_table))]
