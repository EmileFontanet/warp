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
