from scipy.ndimage import gaussian_filter1d
from scipy.signal import medfilt
import numpy as np
import lightkurve as lk
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
from lightkurve.correctors import RegressionCorrector
from lightkurve.correctors import DesignMatrix
from multiprocessing import Process, Queue
import tempfile
import os
import logging
warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO)


def download_and_extract_lcs(star_name, sector=None, limit=5, cutout_size=(35, 35), timeout=45, verbose=True):

    tpfs = load_tpf_by_name(star_name, sector=sector, limit=limit,
                            cutout_size=cutout_size, timeout=timeout, verbose=verbose)
    if tpfs is None:
        return None
    lcs = []
    for tpf in tpfs:
        corrected, uncorrected = lc_from_tpf(tpf, verbose=verbose)
        if corrected is not None and uncorrected is not None:
            lcs.append((corrected, uncorrected))
    return lcs


def lc_from_tpf(tpf, verbose=True):
    try:
        # tpf = tpf.remove_nans()
        finite_mask = np.isfinite(tpf.flux).all(axis=(1, 2))
        tpf = tpf[finite_mask]
        aper = tpf.create_threshold_mask()
        uncorrected_lc = tpf.to_lightcurve(aperture_mask=aper)
        # Make a design matrix and pass it to a linear regression corrector
        dm = DesignMatrix(tpf.flux[:, ~aper], name='regressors').pca(
            5).append_constant()
        rc = RegressionCorrector(uncorrected_lc)
        corrected_ffi_lc = rc.correct(dm)

        # Optional: Remove the scattered light, allowing for the large offset from scattered light
        corrected_ffi_lc = uncorrected_lc - rc.model_lc + \
            np.percentile(rc.model_lc.flux, 5)
    except Exception as e:
        if verbose:
            logging.warning(f"[WARN] Failed to process TPF: {e}")
        return None, None
    return corrected_ffi_lc, uncorrected_lc


def _download_to_fits(row, cutout_size, filename, queue):
    try:
        tpf = row.download(cutout_size=cutout_size)
        tpf.to_fits(filename, overwrite=True)
        queue.put(filename)
    except Exception:
        queue.put(None)


def download_with_timeout(row, cutout_size=(25, 25), timeout=60):
    q = Queue()

    with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as tmp:
        fname = tmp.name

    p = Process(target=_download_to_fits, args=(row, cutout_size, fname, q))
    p.start()
    p.join(timeout)

    if p.is_alive():
        p.terminate()
        p.join()
        return None

    result = q.get() if not q.empty() else None

    if result is None:
        if os.path.exists(fname):
            os.remove(fname)
        return None

    return fname


def load_tpf_by_name(
    star_name,
    sector=None,
    limit=5,
    cutout_size=(25, 25),
    timeout=60,
    verbose=True
):

    if verbose:
        logging.info(f"[INFO] Querying TESS for {star_name}...")

    search_result = lk.search_tesscut(star_name, sector=sector)

    if len(search_result) == 0:
        return None

    if limit is not None:
        search_result = search_result[:limit]

    results = []

    for row in search_result:
        if verbose:
            logging.info(f"Downloading {row.mission}")
            # print(f"Downloading {row.mission}")

        # fname = download_with_timeout(
        #     row, cutout_size=cutout_size, timeout=timeout)

        # if fname is None:
        #     if verbose:
        #         print(f"  [WARN] Timeout for {row.mission}")
        #     continue

        try:
            # tpf = lk.read(fname)
            tpf = row.download(cutout_size=cutout_size)
            results.append(tpf)
            # corrected, uncorrected = lc_from_tpf(tpf)
            # results.append((corrected, uncorrected))
        except Exception as e:
            logging.error(f'Error with {row}')

    return results if len(results) > 0 else None


def to_fits_light(lc, filename, verbose=True):
    from astropy.table import Table

    data = {
        "time": lc.time.value,
        "flux": lc.flux.value,
    }

    if lc.flux_err is not None:
        data["flux_err"] = lc.flux_err.value

    tbl = Table(data)
    tbl.write(filename, format="fits", overwrite=True)
    if verbose:
        logging.info(f"[INFO] Light curve saved to {filename}.")
    return


def detrend(time, flux, flux_err=None,
            smooth_days=3.0,
            sigma_clip=5.0,
            return_trend=False):
    """
    Mild detrending for TESS red-giant light curves.

    Parameters
    ----------
    time : array
        Time in days.
    flux : array
        Flux values.
    flux_err : array or None
        Flux uncertainties.
    smooth_days : float
        Characteristic smoothing timescale in days.
    sigma_clip : float
        Sigma clipping threshold.
    return_trend : bool
        If True, also return fitted trend.

    Returns
    -------
    time_clean, flux_flat[, trend]
    """
    from scipy.interpolate import UnivariateSpline
    from scipy.signal import savgol_filter
    time = np.asarray(time)
    flux = np.asarray(flux)

    mask = np.isfinite(time) & np.isfinite(flux)

    if flux_err is not None:
        flux_err = np.asarray(flux_err)
        mask &= np.isfinite(flux_err)

    t = time[mask]
    f = flux[mask]

    if flux_err is not None:
        ferr = flux_err[mask]
    else:
        ferr = np.ones_like(f)

    # median normalize
    med = np.nanmedian(f)
    f = f / med

    # sigma clipping
    resid = f - np.nanmedian(f)
    mad = 1.4826 * np.nanmedian(np.abs(resid))
    good = np.abs(resid) < sigma_clip * mad

    t = t[good]
    f = f[good]
    ferr = ferr[good]

    # cadence
    dt = np.nanmedian(np.diff(t))

    # smoothing window in cadences
    win = int(np.round(smooth_days / dt))

    if win % 2 == 0:
        win += 1

    win = max(win, 11)

    trend = savgol_filter(f, window_length=win, polyorder=1)

    flux_flat = f / trend - 1.0

    if return_trend:
        return t, flux_flat, trend

    return t, flux_flat


def ls_periodogram(time, flux_flat,
                   fmin=5.0,
                   fmax=150.0,
                   oversample=5):
    """
    Compute Lomb-Scargle periodogram.
    """
    from astropy.timeseries import LombScargle
    fmin_hz = fmin * 1e-6 * 86400.0
    fmax_hz = fmax * 1e-6 * 86400.0

    ls = LombScargle(time, flux_flat)

    freq, power = ls.autopower(
        minimum_frequency=fmin_hz,
        maximum_frequency=fmax_hz,
        samples_per_peak=oversample
    )

    freq_uHz = freq / 86400.0 * 1e6

    return freq_uHz, power


def estimate_numax(freq, power,
                   search=(60, 180),
                   smooth_uHz=8.0,
                   bg_width_uHz=40.0):
    """
    First-pass numax estimate from continuum-normalized PSD.
    """

    mask = (
        np.isfinite(freq)
        & np.isfinite(power)
        & (freq >= search[0])
        & (freq <= search[1])
    )

    f = freq[mask]
    p = power[mask]

    df = np.median(np.diff(f))

    # smooth enough to reveal envelope
    sigma_pix = smooth_uHz / df
    p_smooth = gaussian_filter1d(p, sigma_pix)

    # broad continuum
    k = int(bg_width_uHz / df)
    if k % 2 == 0:
        k += 1
    k = max(k, 21)

    background = medfilt(p_smooth, kernel_size=k)

    # normalize
    p_flat = p_smooth / background

    numax = f[np.argmax(p_flat)]

    return numax, f, p_smooth, background, p_flat, p


def save_for_pysyd(filename, time, flux):
    """
    Save detrended light curve for pySYD.

    Parameters
    ----------
    filename : str
    time : array
        Time in days
    flux : array
        Detrended flux
    """

    arr = np.column_stack((time, flux))

    np.savetxt(
        filename,
        arr,
        fmt="%.8f %.16f"
    )
# def estimate_numax(freq_uHz, power,
#                    fmin=40.0,
#                    fmax=180.0,
#                    smooth_uHz=4.0):
#     """
#     Estimate numax inside a restricted search window.
#     """

#     mask = (
#         np.isfinite(freq_uHz)
#         & np.isfinite(power)
#         & (freq_uHz >= fmin)
#         & (freq_uHz <= fmax)
#     )

#     f = freq_uHz[mask]
#     p = power[mask]

#     df = np.median(np.diff(f))

#     sigma_pix = max(1, smooth_uHz / df)

#     p_smooth = gaussian_filter1d(p, sigma_pix)

#     numax = f[np.argmax(p_smooth)]

#     return numax, f, p_smooth
# def _download_one(row, cutout_size):
#     return row.download(cutout_size=cutout_size)


# def download_tpf_with_timeout(sector, timeout=60, cutout_size=(50, 50), verbose=True):
#     from multiprocessing import Pool

#     with Pool(1) as p:
#         result = p.apply_async(_download_one, (sector, cutout_size))
#         try:
#             return result.get(timeout=timeout)
#         except Exception as e:
#             print(e)
#             print(f'[WARN] TESS download timed out for {sector.mission[0]}.')
#             p.terminate()
#             return None


# def load_tpf_by_name(star_name, sector=None, limit=None, verbose=True, timeout=60):
#     if verbose:
#         print(f"[INFO] Querying TESS for {star_name}...")
#     search_result = lk.search_tesscut(
#         star_name, sector=sector)
#     if len(search_result) == 0:
#         if verbose:
#             print(f"[WARN] No TESS data found for {star_name}.")
#         return None
#     if limit is not None:
#         search_result = search_result[:limit]
#     tpfs = []
#     for sector in search_result:
#         try:
#             if verbose:
#                 print(
#                     f"[INFO] Downloading {sector.mission[0]}...")
#             tpf = download_tpf_with_timeout(
#                 sector, timeout=timeout, cutout_size=(50, 50), verbose=verbose)
#             print(tpf)
#             if tpf is not None:
#                 tpfs.append(tpf)
#         except Exception as e:
#             if verbose:
#                 print(
#                     f"[ERROR] Failed to download TESS data for {star_name}: {e}")

#     return tpfs
