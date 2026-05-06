import requests as rq
from bs4 import BeautifulSoup as bs
import pandas as pd


def get_hip_id(star_name):
    from astroquery.simbad import Simbad

    print(f"[INFO] Retrieving HIP ID for star: {star_name}")
    result_table = Simbad.query_objectids(star_name)
    hip_name = [id for id in result_table['id'] if id.startswith('HIP')][0]
    if len(hip_name) == 0:
        print('[WARN] No HIP name found for this star.')
        return None
    return hip_name


def query_hip_photometry(hip_number):

    url = f'https://cdsarc.cds.unistra.fr/viz-bin/nph-Plot/Vgraph/htm?I/239/{hip_number}&0'
    res = rq.get(url)
    soup = bs(res.text, features="html5lib")
    pre = soup.find('pre')
    a = pre.find('a')
    data = []
    added_header = False
    for line in a.next_sibling.split('\n'):
        if (len(line) == 0):
            continue
        if (added_header == False):
            header = line.split('|')
            header = [c.replace(' ', '') for c in header]
            added_header = True
        else:
            row = {}
            split_line = line.split('|')
            for it, param in enumerate(header):
                row[param] = float(split_line[it])
            data.append(row)
    data = pd.DataFrame(data)
    return data


def query_kervella_table(hip_id):
    from astroquery.vizier import Vizier
    v = Vizier(
        columns=["*"],     # retrieve all columns
        row_limit=50
    )

    result = v.query_constraints(
        catalog="J/A+A/657/A7/tablea1",
        HIP=str(hip_id)
    )
    if len(result) == 0:
        return None
    return result[0].to_pandas()


def query_gaia_nss(gaia_id):
    from astroquery.gaia import Gaia
    query = f"""
     SELECT nss_two_body_orbit.source_id,nss_two_body_orbit.nss_solution_type,nss_two_body_orbit.ra,nss_two_body_orbit.dec,nss_two_body_orbit.parallax,nss_two_body_orbit.pmra,nss_two_body_orbit.pmdec,nss_two_body_orbit.period,nss_two_body_orbit.t_periastron,nss_two_body_orbit.eccentricity,nss_two_body_orbit.center_of_mass_velocity,nss_two_body_orbit.semi_amplitude_primary,nss_two_body_orbit.mass_ratio,nss_two_body_orbit.inclination,nss_two_body_orbit.arg_periastron FROM gaiadr3.nss_two_body_orbit WHERE source_id = {gaia_id}"""

    job = Gaia.launch_job(query)
    return job.get_results().to_pandas()
