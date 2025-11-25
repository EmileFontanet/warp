import numpy as np
import pandas as pd


def mad_clip_mask(values, groups=None, threshold=5, n_iter=3, verbose=False):
    """
    Return boolean mask marking which points are kept by MAD clipping.

    Args:
        values (1D array): Quantity to clip (e.g., RV, FWHM, BIS).
        groups (1D array or None): Optional grouping (e.g., instrument names).
        threshold (float): MAD threshold.
        n_iter (int): Iterations.

    Returns:
        mask (np.array of bool): True = kept, False = clipped.
    """
    values = np.asarray(values)
    mask = np.ones(len(values), dtype=bool)

    if groups is None:
        groups = np.array(["__all__"] * len(values))
    else:
        groups = np.asarray(groups)

    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        submask = np.ones(len(idx), dtype=bool)
        subvals = values[idx]
        if verbose:
            print(f"[INFO] MAD clipping for group: {g}")
        for _ in range(n_iter):
            median = np.median(subvals[submask])
            mad = np.median(np.abs(subvals[submask] - median))
            if mad == 0:
                break

            modz = 0.6745 * (subvals - median) / mad
            new_submask = np.abs(modz) < threshold
            rejected = np.sum(~new_submask & submask)
            if verbose:
                print(f" Iter {_+1}: rejected {rejected} points.")
            # if no change → stop
            if np.all(new_submask == submask):
                break

            submask = new_submask

        mask[idx] = submask

    return mask


def weighted_mean(values, errors):
    """
    Compute weighted mean and its uncertainty.

    Args:
        values (1D array): Values.
        errors (1D array): Uncertainties.
    Returns:
        mean (float): Weighted mean.
        mean_error (float): Uncertainty of the weighted mean.
    """
    weights = 1 / (errors ** 2)
    mean = np.sum(values * weights) / np.sum(weights)
    mean_error = np.sqrt(1 / np.sum(weights))
    return mean, mean_error


def bin_data():
    pass
