from astropy.timeseries import LombScargle
import numpy as np
from .stats import weighted_mean
import pandas as pd


def gls_periodogram(t, y, yerr=None, min_freq=None, max_freq=None, samples=10000):

    ls = LombScargle(t, y, yerr, center_data=True)
    if min_freq is None:
        min_freq = 1 / (max(t) - min(t))
    if max_freq is None:
        dt = np.median(np.diff(np.sort(t)))
        max_freq = 0.5 / dt
        max_freq = max(max_freq, 1.0)
    freq = np.linspace(min_freq, max_freq, samples)
    power = ls.power(freq)

    best_freq = freq[np.argmax(power)]
    best_period = 1 / best_freq
    fap = ls.false_alarm_probability(power.max())
    return freq, power, best_period, fap


def bin_by_night(rv_data, group_cols=['date_night', 'ins_name'], exclude_cols=None, verbose=True):
    if isinstance(group_cols, str):
        group_cols = [group_cols]
    kept_cols = [
        c for c in rv_data.columns if exclude_cols is None or not any(col in c for col in exclude_cols)]
    kept_cols = rv_data[kept_cols].select_dtypes(
        include=[np.number]).columns.tolist()
    if verbose:
        print(f"[INFO] Grouping by columns: {group_cols}")
        print(f"[INFO] Binning columns: {kept_cols}")
        print(
            f"[INFO] Excluding: {[c for c in rv_data.columns if c not in kept_cols]}")
    err_map = {
        col: f"{col}_err"
        for col in kept_cols
        if f"{col}_err" in rv_data.columns
    }
    # We fill the _err entry when we compute the weighted mean of the value
    kept_cols = [c for c in kept_cols if '_err' not in c]
    binned_data = []
    for group_key, group in rv_data.groupby(group_cols):
        if verbose:
            print(f"[INFO] Binning group {group_key}: {len(group)} points")
        row = {}
        # First fill in group key columns
        if isinstance(group_key, tuple):
            for col, val in zip(group_cols, group_key):
                row[col] = val
        else:
            row[group_cols[0]] = group_key

        for col in kept_cols:
            if col in err_map:
                y = group[col].values
                yerr = group[err_map[col]].values
                wmean, wmean_err = weighted_mean(y, yerr)
                row[col] = wmean
                row[err_map[col]] = wmean_err
            elif col == "obj_date_bjd":
                if "spectro_ccf_rv_err" in group:
                    row[col] = weighted_mean(
                        group[col].values,
                        group["spectro_ccf_rv_err"].values
                    )[0]
                else:
                    row[col] = group[col].mean()
            else:
                row[col] = group[col].mean()
        binned_data.append(row)
    binned_df = pd.DataFrame(binned_data)
    return binned_df
