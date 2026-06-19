from .nusyd import nuSYD
import numpy as np
from astropy.io import fits
from .tess import detrend, get_numax_init
import traceback


def process_lc_file(lc_path,  guess_table=None, plot=False, min_nu=None, max_nu=None):
    try:
        f = fits.open(lc_path)
        hdr = f[0].header
        dat = f[1].data
        star = lc_path.split("/")[-2]  # assumes .../HDxxxx/file.fits
        sector = hdr.get("SECTOR", None)
        fname = lc_path.lower()
        if ('QLP' in hdr.get("ORIGIN", "")):
            provenance = "QLP"
        elif ('spoc' in hdr.get("PROCVER", "").lower()):
            provenance = "SPOC"
        else:
            if "qlp" in fname:
                provenance = "QLP"
            elif "spoc" in fname:
                provenance = "SPOC"
            else:
                provenance = "unknown"
        # -------------------------
        # flux selection
        # -------------------------
        if "PDCSAP_FLUX" in dat.columns.names:
            flux_col = "PDCSAP_FLUX"
        elif "KSPSAP_FLUX" in dat.columns.names:
            flux_col = "KSPSAP_FLUX"
        elif "SAP_FLUX" in dat.columns.names:
            flux_col = "SAP_FLUX"
        else:
            flux_col = None
            raise ValueError("No usable flux column found")
        quality_mask = dat['QUALITY'] == 0
        dat = dat[quality_mask]
        time = dat["TIME"]
        flux = dat[flux_col]
        flux_err = dat[flux_col + "_ERR"] if flux_col + \
            "_ERR" in dat.columns.names else None

        time, flux, trend_full = detrend(
            time, flux, flux_err, return_trend=True, smooth_days=2)
        # flux = (flux - np.nanmedian(flux)) / np.nanmedian(flux)

        # -------------------------
        # if guess_table == 'custom':
        #     guess_nu, freq, power, smooth = get_numax_init(time*24*60*60, flux)
        if guess_table is not None:
            guess_nu = guess_table[guess_table.HD ==
                                   star].numax_guess.values[0] if star in guess_table["HD"].tolist() else "from_lc"
        else:
            guess_nu = "from_lc"
        # run nuSYD
        # -------------------------
        runner = nuSYD(
            time,
            flux,
            mc_iter=200,
            lc_type="TESS",
            guess_numax=guess_nu,
            plot=plot
        )

        results = runner.run()
        f.close()

        return {
            "file": lc_path,
            "star": star,
            "sector": sector,
            "provenance": provenance,
            "numax": results.get("numax", np.nan),
            "numax_err": results.get("errors", np.nan),
            'flux': flux_col
        }

    except Exception as e:
        print(e)
        traceback.print_exc()
        return {
            "file": lc_path,
            "error": str(e),
            'star': star,
            'sector': sector,
            'provenance': provenance,
            'flux': flux_col
        }
