from . import cascades
from . import dace
from .simbad import query_simbad
from .hipparcos import get_hip_id, query_hip_photometry
from .stats import mad_clip_mask, weighted_mean
from .tseries import gls_periodogram
import pandas as pd


class Star:
    def __init__(self, name=None, instrument=None, load_rv=True, load_tess=False, max_erv=30,
                 do_adjust_means=True, do_secular_corr=True, skip_ndrs=True, keep_bad_qc=False, verbose=True):
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
            print('loading rv data...')
            self.load_rv()
            if self.do_adjust_means:
                self.adjust_means(verbose=verbose)
        self.ids = {}
        self.coords = None
        self.gaia = None
        self.lc = None
        self.planets = None
        return

    def load_ids(self):
        return

    def load_gaia(self):
        pass

    def load_rv(self):
        self.rv_data = dace.download_points(
            self,
            instrument=self.instrument,
            do_secular_corr=self.do_secular_corr,
            skip_ndrs=self.skip_ndrs
        )
        if not self.keep_bad_qc:
            self.remove_condition(self.rv_data.spectro_drs_qc == False,
                                  origin='drs_qc')
        self.remove_condition(self.rv_data.spectro_ccf_rv_err > self.max_erv,
                              origin=f"rv_err_gt_{self.max_erv}")

    def compute_periodogram(self, min_period=None, max_period=None):
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

    def load_lc(self):
        pass

    def load_planets(self):
        pass

    def load_hip_photometry(self):
        if not hasattr(self, 'hip'):
            self.hip = get_hip_id(self.name)
        if self.hip is None:
            print("[WARN] No HIP ID found, cannot retrieve photometry.")
            return None
        self.hip_photometry = query_hip_photometry(hip_number=self.hip[4:])
        return self.hip_photometry

    def clip_rv(self, threshold=5, n_iter=3, verbose=False):
        if not hasattr(self, 'rv_data'):
            print("[WARN] No RV data loaded, cannot perform MAD clipping.")
            return
        mask = mad_clip_mask(
            self.rv, groups=self.rv_data.ins_name, threshold=threshold, n_iter=n_iter,
            verbose=verbose)
        self.remove_condition(~mask, origin=f'rv_mad_clip_{threshold}sigma')
        # self.rv_data = self.rv_data[mask].copy()
        return

    def adjust_means(self, verbose=True):
        if not hasattr(self, 'rv_data'):
            print("[WARN] No RV data loaded, cannot adjust means.")
            return
        for ins in self.rv_data.ins_name.unique():
            ins_mask = self.rv_data.ins_name == ins
            mean_rv, _ = weighted_mean(self.rv[ins_mask],
                                       self.rv_err[ins_mask])
            self.rv_data.loc[ins_mask, 'spectro_ccf_rv'] -= mean_rv
            if verbose:
                print(
                    f"[INFO] Adjusted mean RV for instrument {ins} by subtracting {mean_rv:.3f} m/s.")
        self.did_adjust_means = True
        return

    def remove_condition(self, condition, origin, verbose=True):
        """
        Remove rows from rv_data based on a boolean condition.

        Parameters
        ----------
        condition : array-like of bool
            True for rows to remove.
        origin : str
            Reason string to store in bad_points['reason'].
        """

        # ensure boolean mask length matches
        if len(condition) != len(self.rv_data):
            raise ValueError(
                "Condition mask length does not match rv_data length.")

        # extract rows to remove
        removed = self.rv_data[condition].copy()
        removed["reason"] = origin
        if verbose:
            n_removed = len(removed)
            print(
                f"[INFO] Removing {n_removed} points from rv_data for reason: {origin}")
        # append to bad_points
        self.bad_points = pd.concat(
            [self.bad_points, removed], ignore_index=True)

        # keep only remaining rows
        self.rv_data = self.rv_data[~condition].reset_index(drop=True)

    @property
    def rv(self):
        return self.rv_data.spectro_ccf_rv if self.rv_data is not None else None

    @property
    def rv_err(self):
        return self.rv_data.spectro_ccf_rv_err if self.rv_data is not None else None

    @property
    def t(self):
        return self.rv_data.obj_date_bjd if self.rv_data is not None else None
