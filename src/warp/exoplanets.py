import requests as rq
import pandas as pd


def fetch_nasa_archive(table='ps', only_default=True, columns=None, params=None):
    """
    Fetch exoplanet data from NASA Exoplanet Archive.

    Args:
        table (str): Table name to query (default: 'ps' for Planetary Systems).
        only_default (bool): Whether to fetch only default parameters for each systems.
        columns (list or str or None): Specific columns to retrieve. Default is None, which retrieves only planet names. If 'main', retrieves a predefined set of main parameters.
        params (dict or None): Additional query parameters.

    Returns:
        pd.DataFrame: DataFrame containing the exoplanet data.
    """

    base_url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query='
    if columns is not None:
        if isinstance(columns, list):
            columns = ','.join(columns)
        elif isinstance(columns, str):
            if (columns == 'main'):
                columns = 'pl_name,hostname,pl_orbper,pl_controv_flag,discoverymethod,pl_refname,pl_massj, pl_bmassj,pl_msinij,st_refname,st_spectype,st_teff,st_rad,st_mass,st_logg'
            columns = columns
        query = f"SELECT+{columns}+FROM+{table}+where"

    else:
        query = f"SELECT+pl_name+FROM+{table}+where"

    if only_default:
        query += "+default_flag=1&"
    else:
        query += "+1=1&"
    if params is not None:
        for key, value in params.items():
            query += f"+{key}='{value}'+&"
    query = query+'format=json'
    response = rq.get(base_url+query)
    response.raise_for_status()
    df = pd.DataFrame(response.json())

    return df


def get_columns(table='ps'):
    """
    Get available columns for a given table in the NASA Exoplanet Archive.

    Args:
        table (str): Table name (default: 'ps' for Planetary Systems).

    Returns:
        list: List of column names.
    """
    query = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query= select+table_name,column_name,description,datatype,column_index+from+TAP_SCHEMA.columns+%20where+table_name+like+ %27ps%27&format=json"
    response = rq.get(query)
    response.raise_for_status()
    df = pd.DataFrame(response.json())
    return df


def is_planet_published(pl_names, planet_list=None):
    from astroquery.simbad import Simbad

    """
    Check if a planet is published in the NASA Exoplanet Archive.

    Args:
        pl_names (str or list): Name(s) of the planet(s).
        planet_list (pd.DataFrame or None): Optional pre-fetched DataFrame of planets to check against. If None, the function will fetch the data from the archive.
    Returns:
        list containing the pl parameters and refs if the planet is published, False otherwise.
        """
    Simbad.reset_votable_fields()
    Simbad.add_votable_fields('ids')
    if isinstance(pl_names, str):
        pl_names = [pl_names]
    try:
        result_table = Simbad.query_objects(pl_names).to_pandas()
    except Exception as e:
        print('Could not retrieve SIMBAD IDs:', e)
        return None
    if planet_list is None:
        all_planets = fetch_nasa_archive(
            columns=['pl_name', 'pl_orbper', 'pl_refname', 'hostname'], only_default=False)
    else:
        all_planets = planet_list
    results = []
    for it, planet_name in enumerate(pl_names):
        temp_result = []
        ids = result_table.iloc[it]['ids']
        ids_low = [c.lower().replace(" ", "") for c in ids.split('|')]
        matches = all_planets[all_planets.hostname.str.lower(
        ).str.replace(" ", "").isin(ids_low)]
        if len(matches) > 0:
            for i in range(len(matches)):
                temp_result.append({
                    'pl_name': matches.iloc[i]['pl_name'],
                    'pl_orbper': matches.iloc[i]['pl_orbper'],
                    'pl_refname': matches.iloc[i]['pl_refname']
                })
            results.append(temp_result)
        else:
            results.append(None)
    # planet_data = all_planets[all_planets['pl_name'] == pl_name]
    # if len(planet_data) == 0:
    #     return False
    return results
