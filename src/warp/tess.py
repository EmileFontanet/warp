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


def download_and_extract_lcs(star_name, sector=None, limit=35, cutout_size=(35, 35), timeout=45, verbose=True, download='individual'):
    """Downloads and extracts light curves of the TESS mission from the tesscut images

    Args:
        star_name (str): Name of the star to download the light curve for.
        sector (list, optional): List of sectors to download. Defaults to None.
        limit (int, optional): Maximum number of TPFs to download. Defaults to 35.
        cutout_size (tuple, optional): Size of the cutout image. Defaults to (35, 35).
        timeout (int, optional): Timeout for the download. Defaults to 45.
        verbose (bool, optional): Whether to print verbose output. Defaults to True.
        download (str, optional): Type of download to perform. Defaults to 'individual'. Options are 'individual' or 'bulk'.

    Returns:
        _type_: _description_
    """
    tpfs = load_tpf_by_name(star_name, sector=sector, limit=limit,
                            cutout_size=cutout_size, timeout=timeout, verbose=verbose, download=download)
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
        uncorrected_lc = uncorrected_lc.remove_nans()
        # Make a design matrix and pass it to a linear regression corrector
        dm = DesignMatrix(tpf.flux[:, ~aper], name='regressors').pca(
            5).append_constant()
        rc = RegressionCorrector(uncorrected_lc)
        corrected_ffi_lc = rc.correct(dm)

        # Optional: Remove the scattered light, allowing for the large offset from scattered light
        corrected_ffi_lc = uncorrected_lc - rc.model_lc + \
            np.percentile(rc.model_lc.flux, 5)
        corrected_ffi_lc = corrected_ffi_lc.remove_nans()
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
    limit=35,
    cutout_size=(25, 25),
    timeout=60,
    verbose=True,
    download='bulk'

):

    if verbose:
        logging.info(f"[INFO] Querying TESS for {star_name}...")

    search_result = lk.search_tesscut(star_name, sector=sector)

    if len(search_result) == 0:
        return None

    if limit is not None:
        search_result = search_result[:limit]

    results = []
    if download == 'bulk':
        try:
            results = search_result.download_all()
            return results if len(results) > 0 else None
        except Exception as e:
            logging.error(f'Error occurred while downloading all TPFs: {e}')
            return None
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
            return_trend=False,
            verbose=True,
            goal='astero'):
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
    goal : str
        The goal of the detrending. Options are 'astero' or 'transit'.

    Returns
    -------
    time_clean, flux_flat[, trend]
    """
    from scipy.interpolate import UnivariateSpline, interp1d
    from scipy.signal import savgol_filter
    time = np.asarray(time)
    flux = np.asarray(flux)
    if goal not in ['astero', 'transit']:
        raise ValueError("Invalid goal. Must be 'astero' or 'transit'.")
    mask = np.isfinite(time) & np.isfinite(flux) & (flux > 0)

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
    if (goal == 'astero'):
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

        flux_flat = f / trend

        if return_trend:
            return t, flux_flat, trend

        return t, flux_flat
    elif (goal == 'transit'):
        trend_good = np.abs(resid) < 3 * mad
        lc_good = (resid < 10 * mad)  # & (f > 0.93) & (f < 1.05)

        t_lc = t[lc_good]
        f_lc = f[lc_good]
        ferr_lc = ferr[lc_good]

        t_trend = t[trend_good]
        f_trend = f[trend_good]

        # cadence
        dt = np.nanmedian(np.diff(t_trend))

        # smoothing window in cadences
        win = int(np.round(smooth_days / dt))

        if win % 2 == 0:
            win += 1

        win = max(win, 11)

        trend = savgol_filter(f_trend, window_length=win, polyorder=1)
        interp = interp1d(t_trend, trend,
                          bounds_error=False,
                          fill_value='extrapolate')

        trend_full = interp(t_lc)
        flux_flat = f_lc / trend_full

        if return_trend:
            return t_lc, flux_flat, trend_full

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


def get_predicted_transits(rvmod, lightcurve):
    """Get the predicted transit times based on the RV fit

    Args:
        rvmod (rv_model): The rv_model
        lightcurve (hdul): The hdul containing the lightcurve data

    Returns:
        _type_: _description_
    """
    for i in range(rvmod.nkep):
        rvmod.set_keplerian_param(f'{i}', ['P', 'e', 'w', 'K', 'Tc'])

    future_transit_times = {}
    future_transit_err = {}
    lc_data = lightcurve[1].data
    for pla in range(rvmod.nkep):
        planet_idx = pla
        period_err = rvmod.get_param_error(
        )[1][rvmod.fit_param.index(f'kep.{planet_idx}.P')]
        Tc_err = rvmod.get_param_error(
        )[1][rvmod.fit_param.index(f'kep.{planet_idx}.Tc')]
        transit_time = rvmod.get_param(
            f'kep.{planet_idx}.Tc') + rvmod.t0 + 2400000
        future_transit_times[f'planet_{pla+1}'] = []
        future_transit_err[f'planet_{pla+1}'] = []
        lc_times = lc_data['time']+2457000
        n_start = int((lc_times[0]-transit_time) /
                      rvmod.get_param(f'kep.{planet_idx}.P'))-1
        # future_transit_times[f'planet_{pla+1}'][f'sector_{sect}'] = []
        # future_transit_err[f'planet_{pla+1}'][f'sector_{sect}'] = []
        for n in np.linspace(n_start, n_start+499, 500, dtype='int'):
            future_transit = transit_time + \
                (n*rvmod.get_param(f'kep.{planet_idx}.P'))
            transit_err = np.sqrt((n*period_err)**2 + Tc_err**2)
            future_transit_lower = future_transit-transit_err
            future_transit_upper = future_transit+transit_err
            if (future_transit > lc_times[0] and future_transit < lc_times[-1]) or (future_transit_lower > lc_times[0] and future_transit_lower < lc_times[-1]) or (future_transit_upper > lc_times[0] and future_transit_upper < lc_times[-1]):
                future_transit_times[f'planet_{pla+1}'].append(
                    future_transit)
                future_transit_err[f'planet_{pla+1}'].append(
                    (future_transit_lower, future_transit_upper))

    return future_transit_times, future_transit_err


def list_tess_lightcurves(target, output_format='pandas', accepted_pipelines=['SPOC', 'TESS-SPOC', 'QLP'], verbose=True):
    from astroquery.mast import Observations
    """
    target: TIC ID string ('TIC 123456789') or coordinates
    output_format: str, either 'pandas' or 'original'
    accepted_pipelines: list of str, the pipelines to accept (by order of preference)
    """

    obs = Observations.query_criteria(
        target_name=target,
        dataproduct_type="timeseries"
    )

    # Filter by accepted pipelines
    mask = obs['provenance_name'] == None
    for pipeline in accepted_pipelines:
        mask |= obs['provenance_name'] == pipeline
    obs = obs[mask]

    if len(obs) == 0:
        print(f"No light curves found for {target}")
        return None
    for sector in np.unique(obs['sequence_number']):
        if len(obs[obs['sequence_number'] == sector]) > 1:
            if verbose:
                logging.warning(
                    f"Multiple light curves found for {target} in sector {sector}. Attempting to filter by pipeline.")
            for pipeline in accepted_pipelines:
                if pipeline in obs['provenance_name'][obs['sequence_number'] == sector]:
                    obs = obs[~((obs['sequence_number'] == sector)
                                & (obs['provenance_name'] != pipeline))]
                    break

    if output_format == 'pandas':
        df = obs.to_pandas()
        cols = [
            "target_name",
            "sequence_number",
            "provenance_name",
            "obs_id",
            "t_exptime",
            "obs_title"
        ]

        df = df[cols].sort_values(
            ["sequence_number", "provenance_name"]
        )
    elif output_format == 'both':
        df = obs.to_pandas()
        cols = [
            "target_name",
            "sequence_number",
            "provenance_name",
            "obs_id",
            "t_exptime",
            "obs_title"
        ]

        df = df[cols].sort_values(
            ["sequence_number", "provenance_name"]
        )
        return df, obs
    else:
        return obs

    return df


def get_numax_init(time, flux):
    from scipy.signal import savgol_filter
    from astropy.timeseries import LombScargle
    # remove NaNs
    m = np.isfinite(time) & np.isfinite(flux)
    time = time[m]
    flux = flux[m]

    flux = flux - np.nanmedian(flux)
    flux = flux / np.nanmedian(np.abs(flux))
    # -------------------------
    # PSD (Lomb-Scargle)
    # -------------------------
    freq = np.linspace(1, 300, 5000)  # µHz region of interest

    ls = LombScargle(time, flux)
    power = ls.power(freq * 1e-6)  # convert µHz → Hz internally

    power = np.array(power)

    # -------------------------
    # -------------------------
    smooth = savgol_filter(power, 151, 3)

    # -------------------------
    # restrict to red giant region
    # -------------------------
    mask = (freq >= 10) & (freq <= 250)

    f_sub = freq[mask]
    p_sub = smooth[mask]
    # νmax guess = peak
    nu_max = f_sub[np.argmax(p_sub)]
    print(nu_max)
    return nu_max, freq, power, smooth
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
