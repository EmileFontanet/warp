from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy import units as u


def query_gaia(name, radius=2.0*u.arcsec):
    coord = SkyCoord(name) if isinstance(
        name, str) and ',' in name else SkyCoord.from_name(name)
    j = Gaia.cone_search_async(coord, radius=radius).get_results()
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
