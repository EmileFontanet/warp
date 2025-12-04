import requests as rq
import pandas as pd


def fetch_nasa_archive(table='ps', only_default=True, columns=['pl_name'], params=None):
    """
    Fetch exoplanet data from NASA Exoplanet Archive.

    Args:
        table (str): Table name to query (default: 'ps' for Planetary Systems).
        only_default (bool): Whether to fetch only default parameters for each systems.
        columns (list or str or None): Specific columns to retrieve.
        params (dict or None): Additional query parameters.

    Returns:
        pd.DataFrame: DataFrame containing the exoplanet data.
    """
    base_url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query='
    if columns is not None:
        if isinstance(columns, list):
            columns = ','.join(columns)
        elif isinstance(columns, str):
            columns = columns
        query = f"SELECT+{columns}+FROM+{table}+where"

    else:
        query = f"SELECT+*+FROM+{table}+where"

    if only_default:
        query += "+default_flag=1&"

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
