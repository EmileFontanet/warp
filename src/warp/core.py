from . import cascades
from . import dace
from .simbad import query_simbad
from .hipparcos import get_hip_id, query_hip_photometry
from .stats import mad_clip_mask, weighted_mean
from .tseries import gls_periodogram
import pandas as pd
from .config import accepted_pipelines
import numpy as np


class Star:
    def __init__(self, name=None, instrument=None, load_rv=True, load_tess=False, max_erv=30,
                 do_adjust_means=True, do_secular_corr=True, skip_ndrs=True, keep_bad_qc=False, verbose=True,
                 filter_pipeline=True, filter_columns=True, remove_negative_erv=True):
        """A class instance representing a star. This object allows to download, visualize and analyse RV and 
        photometric data.

        Args:
            name (str, optional): The name of the target. Defaults to None.
            instrument (str or list, optional): The instruments for which the RV data should be downloaded. Defaults to None.
            load_rv (bool, optional): Whether to download the RV data from DACE or not. Defaults to True.
            load_tess (bool, optional): Wheter to download data from TESS. Defaults to False.
            max_erv (float, optional): The maximum RV uncertainty accepted. Points with a larger error will be discarded. Defaults to 30.
            do_adjust_means (bool, optional): If true, the weighted average per instrument will be subtracted from the data. Defaults to True.
            do_secular_corr (bool, optional): If true, the data will be corrected for secular acceleration using Gaia or Simbad. Defaults to True.
            skip_ndrs (bool, optional): If true, points obtained with the new CORALIE drs will be discarded. Defaults to True.
            keep_bad_qc (bool, optional): If true, points which fail the drs QC are kept. Defaults to False.
            verbose (bool, optional): If False, all the printing/loggingis suppressed. Defaults to True.
            filter_pipeline (bool, optional): If true, only the latest known pipeline is kept for each instrument. Defaults to True.
            filter_columns (bool, optional): If true, the columns for the RV data obtained from DACE are filtered in order to keep only relevant data. Defaults to True.
            remove_negative_erv (bool, optional): If true, points with a negative RV uncertainty are discarded. Defaults to True.
        """
        self.name = name
        self.instrument = instrument
        self.max_erv = max_erv
        self.do_adjust_means = do_adjust_means
        self.do_secular_corr = do_secular_corr
        self.skip_ndrs = skip_ndrs
        self.keep_bad_qc = keep_bad_qc
        self.did_adjust_means = False
        self.bad_points = pd.DataFrame()
        if name is not None and load_rv:
            self.load_rv(filter_columns=filter_columns,
                         filter_pipeline=filter_pipeline,
                         remove_negative_erv=remove_negative_erv, verbose=verbose)
            if self.do_adjust_means:
                self.adjust_means(verbose=verbose)

        self._ids = None
        self.coords = None
        self.gaia = None
        self.lc = None
        self.planets = None
        return

    def load_ids(self):
        from .simbad import get_ids
        self._ids = get_ids(self.name)
        return

    def load_rv(self, filter_columns=True,
                filter_pipeline=True,
                remove_negative_erv=True, verbose=True):
        self.rv_data = dace.download_points(
            self,
            instrument=self.instrument,
            do_secular_corr=self.do_secular_corr,
            skip_ndrs=self.skip_ndrs,
            verbose=verbose
        )
        self.rv_data = self.rv_data.sort_values(by='obj_date_bjd')
        if filter_columns:
            self.filter_columns()
        if filter_pipeline:
            self.filter_pipeline(verbose=verbose)
        if not self.keep_bad_qc:
            self.remove_condition(self.rv_data.spectro_drs_qc == False,
                                  origin='drs_qc', verbose=verbose, adjust_means=False)
        self.remove_condition(self.rv_data.spectro_ccf_rv_err > self.max_erv,
                              origin=f"rv_err_gt_{self.max_erv}", verbose=verbose, adjust_means=False)
        if remove_negative_erv:
            self.remove_condition(self.rv_data.spectro_ccf_rv_err <= 0,
                                  origin='negative_rv_err', verbose=verbose, adjust_means=False)

    def filter_pipeline(self, verbose=True):
        from .utils import get_latest_pipeline
        if self.rv_data is None:
            print("[WARN] No RV data loaded, cannot filter pipelines.")
            return
        for ins in self.ins.unique():
            accepted_pipeline = get_latest_pipeline(
                ins, self.rv_data[self.ins == ins]['ins_drs_version'].unique(), verbose=verbose)
            rejected_mask = (self.ins == ins) & (
                ~self.rv_data['ins_drs_version'].isin(accepted_pipeline))
            self.rv_data = self.rv_data[~rejected_mask].copy()
        return

    def filter_columns(self):
        from .config import ignored_cols
        kept_cols = [col for col in self.rv_data.columns if
                     not any(c in col for c in ignored_cols)]
        self.rv_data = self.rv_data[kept_cols].copy()
        return

    def compute_periodogram(self, min_period=None, max_period=None):
        """Computes the periodogram of the RV time series

        Args:
            min_period (float, optional): The minimum period in days. Defaults to None.
            max_period (float, optional): The maximum period in days. Defaults to None.

        Returns:
            Tuple: Returns a tuple (freq, power, best_period, fap) containing the frequency array,
                   power spectrum, best period, and false alarm probability of the highest peak.
        """
        if self.did_adjust_means is False:
            print(
                "[WARN] Means have not been adjusted yet. Consider adjusting means before computing periodogram.")
        if self.rv_data is None:
            print("[WARN] No RV data loaded, cannot compute periodogram.")
            return None
        min_freq = 1 / max_period if max_period is not None else None
        max_freq = 1 / min_period if min_period is not None else None
        freq, power, best_period, fap = gls_periodogram(self.t,
                                                        self.rv,
                                                        self.rv_err,
                                                        min_freq=min_freq,
                                                        max_freq=max_freq)
        return freq, power, best_period, fap

    def load_hip_photometry(self):
        if not hasattr(self, 'hip'):
            self.hip = get_hip_id(self.name)
        if self.hip is None:
            print("[WARN] No HIP ID found, cannot retrieve photometry.")
            return None
        self.hip_photometry = query_hip_photometry(hip_number=self.hip[4:])
        return self.hip_photometry

    def clip_rv(self, threshold=5, n_iter=3, groups='ins_name', ins_list=None, verbose=True, adjust_means=True, inplace=True):
        if not hasattr(self, 'rv_data'):
            print("[WARN] No RV data loaded, cannot perform MAD clipping.")
            return
        mask = mad_clip_mask(
            self.rv, groups=self.rv_data[groups], threshold=threshold, n_iter=n_iter,
            verbose=verbose)
        if (inplace):
            if (~mask).sum() > 0:
                self.remove_condition(
                    ~mask, origin=f'rv_mad_clip_{threshold}sigma', verbose=verbose, adjust_means=adjust_means)

        # self.rv_data = self.rv_data[mask].copy()
        return mask

    def adjust_means(self, verbose=True, ins_list=None, series=None):
        if series is None:
            series = ['spectro_ccf_rv',
                      'spectro_ccf_bispan', 'spectro_ccf_fwhm', 'spectro_ccf_contrast']
        if ins_list is None:
            ins_list = self.rv_data.ins_name.unique()
        if not hasattr(self, 'rv_data'):
            print("[WARN] No RV data loaded, cannot adjust means.")
            return

        for ins in ins_list:
            ins_mask = self.rv_data.ins_name == ins
            mean_rv, _ = weighted_mean(self.rv[ins_mask],
                                       self.rv_err[ins_mask])
            self.rv_data.loc[ins_mask, 'spectro_ccf_rv'] -= mean_rv
            if verbose:
                print(
                    f"[INFO] Adjusted mean RV for instrument {ins} by subtracting {mean_rv:.3f} m/s.")
        self.did_adjust_means = True
        return

    def remove_condition(self, condition, origin, verbose=True, adjust_means=True):
        """
        Remove rows from rv_data based on a boolean condition.

        Parameters
        ----------
        condition : array-like of bool
            True for rows to remove.
        origin : str
            Reason string to store in bad_points['reason'].
        """
        if condition.sum() == 0:
            if verbose:
                print(f"[INFO] No points to remove for reason: {origin}")
            return

        # ensure boolean mask length matches
        if len(condition) != len(self.rv_data):
            raise ValueError(
                "Condition mask length does not match rv_data length.")

        # extract rows to remove
        removed = self.rv_data[condition].copy()
        removed["origin"] = origin
        if verbose:
            n_removed = len(removed)
            print(
                f"[INFO] Removing {n_removed} points from rv_data for reason: {origin}")
        # append to bad_points
        self.bad_points = pd.concat(
            [self.bad_points, removed], ignore_index=True)

        # keep only remaining rows
        self.rv_data = self.rv_data[~condition].reset_index(drop=True)
        if adjust_means:
            self.adjust_means(verbose=verbose)

    def plot_rv(self, ax=None, fig=None, save_fig=False, show_plot=False, **kwargs):
        from .plotting import plot_rv
        fig, ax = plot_rv(self.rv_data, ax=ax, fig=fig, star_name=self.name,
                          save_fig=save_fig, show_plot=show_plot, **kwargs)
        return fig, ax

    def plot_series(self, quantity, ax=None, fig=None, save_fig=False, show_plot=False, **kwargs):
        from .plotting import plot_series
        if quantity not in self.rv_data.columns:
            raise ValueError(
                f"Quantity {quantity} not found in rv_data columns.")
        fig, ax = plot_series(self.rv_data, quantity, ax=ax, fig=fig,
                              star_name=self.name,
                              save_fig=save_fig, show_plot=show_plot, **kwargs)
        return fig, ax

    def bin_by_night(self, group_cols=['date_night', 'ins_name'], exclude_cols=None, verbose=True):
        from .tseries import bin_by_night
        from .config import exclude_cols_bin_by_night
        if exclude_cols is None:
            exclude_cols = exclude_cols_bin_by_night
        self.rv_data = bin_by_night(
            self.rv_data,
            group_cols=group_cols,
            exclude_cols=exclude_cols,
            verbose=verbose
        )
        return

    def fit_keplerian(self, N_pla=3, n_lin=1, stellar_jitter=0, fap_threshold=1e-3, periods_init=[]):
        from .kepmodel_wrapper import fit_keplerian
        self.rv_model = fit_keplerian(
            self.t.values,
            self.rv.values,
            self.rv_err.values,
            self.ins,
            N_pla=N_pla,
            n_lin=n_lin,
            stellar_jitter=stellar_jitter,
            fap_threshold=fap_threshold,
            periods_init=periods_init
        )
        return self.rv_model

    def plot_keplerian_fit(self, ins_list=None, ax=None, fig=None, save_fig=False, show_plot=False, **kwargs):
        from .plotting import plot_rv_keplerian_fit
        if not hasattr(self, 'rv_model'):
            raise ValueError("No rv_model found. Fit a Keplerian model first.")
        if ins_list is not None:
            rv_data = self.rv_data[self.rv_data.ins_name.isin(ins_list)]
        else:
            rv_data = self.rv_data
        fig, ax = plot_rv_keplerian_fit(
            rv_data,
            self.rv_model,
            ax=ax,
            fig=fig,
            star_name=self.name,
            save_fig=save_fig,
            show_plot=show_plot,
            **kwargs
        )
        return fig, ax

    @property
    def rv(self):
        return self.rv_data.spectro_ccf_rv if self.rv_data is not None else None

    @property
    def rv_err(self):
        return self.rv_data.spectro_ccf_rv_err if self.rv_data is not None else None

    @property
    def raw_f(self):
        return self.rv_data.file_rootpath if 'file_rootpath' in self.rv_data.columns else None

    @property
    def t(self):
        return self.rv_data.obj_date_bjd if self.rv_data is not None else None

    @property
    def ins(self):
        return self.rv_data.ins_name if self.rv_data is not None else None

    @property
    def n_points(self):
        return len(self.rv_data) if self.rv_data is not None else 0

    @property
    def n_points_by_ins(self):
        if self.rv_data is None:
            return None
        n_points_dict = {}
        for ins in self.rv_data.ins_name.unique():
            ins_mask = self.rv_data.ins_name == ins
            n_points_dict[ins] = np.sum(ins_mask)
        return n_points_dict

    @property
    def rms(self):
        return np.std(self.rv) if self.rv_data is not None else None

    @property
    def mad(self):
        if self.rv_data is None:
            return None
        median = np.median(self.rv)
        mad = np.median(np.abs(self.rv - median))
        return mad

    @property
    def mad_by_ins(self):
        if self.rv_data is None:
            return None
        mad_dict = {}
        for ins in self.rv_data.ins_name.unique():
            ins_mask = self.rv_data.ins_name == ins
            median = np.median(self.rv[ins_mask])
            mad = np.median(np.abs(self.rv[ins_mask] - median))
            mad_dict[ins] = mad
        return mad_dict

    @property
    def rms_by_ins(self):
        if self.rv_data is None:
            return None
        rms_dict = {}
        for ins in self.rv_data.ins_name.unique():
            ins_mask = self.rv_data.ins_name == ins
            rms_dict[ins] = np.std(self.rv[ins_mask])
        return rms_dict

    @property
    def tspan(self):
        return self.t.max() - self.t.min() if self.rv_data is not None else None

    @property
    def mean(self):
        return np.mean(self.rv) if self.rv_data is not None else None

    @property
    def wmean(self):
        return weighted_mean(self.rv, self.rv_err)[0] if self.rv_data is not None else None

    @property
    def median(self):
        return np.median(self.rv) if self.rv_data is not None else None

    @property
    def ids(self):
        if self._ids is None:
            self.load_ids()
        return self._ids
