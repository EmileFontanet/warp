import pandas as pd
import numpy as np
import os

prob_type1 = [
    'HD8651',
    'HD67929',
    'HD79202',
    'HD112410',
    'HD207229',
    'HD212953',
    'HD58961',
    'HD18448',
    'HD118319',
    'HD224949',
    'HD213770',
    'HD219507',
    'HD199845',
    'HD208285']
prob_type2 = {
    'HD28732': 3/1,
    'HD44880':7/5,
    'HD46099':1.5,
    'HD55964':4/3,
    'HD64121':7/5,
    'HD78271':5/4,
    'HD94890':3/1,
    'HD127195':8/5,
    'HD128853':3/2,
    'HD132905': 3,
       
}
prob_type3= ['HD175219']

def load_star_params():
    """
    Load stellar parameters from a CSV file into a pandas DataFrame.

    Returns:
        pd.DataFrame: DataFrame containing stellar parameters.
    """
    # Assuming the CSV file is named 'star_params.csv' and is in the same directory
    cols = ['HD', 'Stype', 'V', 'V_e', 'B-V', 'B-V_e', 'BC', 'BC_e', 'pi', 'pi_e', 'd', 'd_e', 'd_E', 'Mv', 'Mv_e', 'Bp-Rp', 'Bp-Rp_e', 'G', 'G_e', 'Teff', 'Teff_err', 'logg','logg_e', 'feh', 'feh_e', 'Mstar', 'Mstar_e', 'Lstar', 'Lstar_e', 'Rstar', 'Rstar_e']
    file_path = os.path.join(os.path.dirname(__file__), "data", "cascades.dat")
    df = pd.read_csv(file_path, sep='\s+', on_bad_lines='skip', names=cols)
    df['HD']= np.array(['HD'+str(c) for c in df['HD']])
    return df
