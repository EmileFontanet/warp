import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
ins_colors = {
    'CORALIE98': 'tab:blue',
    'CORALIE07': 'tab:orange',
    'CORALIE14': 'tab:green',
    'CORALIE24': 'tab:red',
    'HARPS03': 'tab:purple',
    'HARPS15': 'tab:brown',
    'NIRPS': 'tab:pink',
    'ESPRESSO19': 'tab:gray',
    'ESPRESSO18': 'tab:olive',
    'CORALIE98-NDRS': 'darkred',
    'CORALIE07-NDRS': 'gold',
    'CORALIE14-NDRS': 'darkblue',
    'CORALIE24-NDRS': 'mediumaquamarine',
    'default': 'tab:cyan'
}


def plot_rv(rv_data, ax=None, fig=None, star_name=None, save_fig=False, show_plot=False, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10, 6)))
    for ins in rv_data['ins_name'].unique():
        for drs in rv_data[rv_data.ins_name == ins]['ins_drs_version'].unique():
            ins_mask = (rv_data['ins_name'] == ins) & (
                rv_data['ins_drs_version'] == drs)
            color = ins_colors.get(ins, ins_colors['default'])
            ax.errorbar(rv_data.loc[ins_mask, 'obj_date_bjd'],
                        rv_data.loc[ins_mask, 'spectro_ccf_rv'],
                        yerr=rv_data.loc[ins_mask, 'spectro_ccf_rv_err'],
                        fmt='o', label=f"{ins} ({drs})", color=color, alpha=0.7,
                        capsize=3, **kwargs)
    ax.set_xlabel('Time [d]')
    ax.set_ylabel('Radial Velocity [m/s]')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_title(
        f'Radial Velocity Measurements for {star_name}' if star_name else "Radial Velocity Measurements")
    if save_fig and fig is not None:
        fig.savefig(f"{star_name}_rv.png", dpi=300)
    if show_plot:
        plt.show(block=False)
    return fig, ax


def plot_rv_keplerian_fit(rv_data, rv_model, ax=None, fig=None, star_name=None,
                          save_fig=False, show_plot=False, date_format='rjd', legend=True, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10, 6)))

    for ins in rv_data['ins_name'].unique():
        ins_mask = rv_data['ins_name'] == ins
        color = ins_colors.get(ins, ins_colors['default'])
        if f"lin.{ins}_offset" in rv_model.fit_param:
            offset = rv_model.get_param(
                f"lin.{ins}_offset")
        elif hasattr(rv_model, 'fixed_offsets') and ins in rv_model.fixed_offsets:
            offset = rv_model.fixed_offsets[ins] + \
                rv_model.get_param('lin.rv_CORALIE_offset')
        elif hasattr(rv_model, 'fixed_offsets') and 'CORALIE14' in ins:
            offset = rv_model.get_param('lin.rv_CORALIE_offset')
        else:
            offset = 0.0
        if date_format == 'datetime':
            time_vals = pd.to_datetime(rv_data.loc[ins_mask, 'date_night'])
            ax.set_xlabel('Time [d]')
        elif date_format == 'rjd':
            time_vals = rv_data.loc[ins_mask, 'obj_date_bjd']
            ax.set_xlabel('Date')
        ax.errorbar(time_vals,
                    rv_data.loc[ins_mask, 'spectro_ccf_rv'] - offset,
                    yerr=rv_data.loc[ins_mask, 'spectro_ccf_rv_err'],
                    fmt='o', label=ins, color=color, alpha=0.7,
                    capsize=2, ms=2, **kwargs)
    # Plot model
    time_grid = np.linspace(np.min(rv_data['obj_date_bjd']),
                            np.max(rv_data['obj_date_bjd']), 5000)
    model_rv = rv_model.keplerian_model(
        time_grid-np.median(rv_data['obj_date_bjd']))
    if date_format == 'datetime':
        time_grid = pd.date_range(start=pd.to_datetime(rv_data['date_night']).min(),
                                  end=pd.to_datetime(
                                      rv_data['date_night']).max(),
                                  periods=5000)
    drift_pow = len([c for c in rv_model.fit_param if 'drift_pow' in c])
    drift = np.zeros_like(time_grid)
    for k in range(drift_pow):
        drift += rv_model.get_param(
            f"lin.drift_pow{k+1}") * (time_grid - np.median(rv_data['obj_date_bjd']))**(k+1)
    ax.plot(time_grid, model_rv + drift,
            color='black', lw=1, label='Keplerian Fit')
    if drift_pow > 0:
        ax.plot(time_grid, drift, color='gray',
                lw=1, ls='--', label=f'Degree {drift_pow} Drift')
    else:
        ax.axhline(0, color='gray', lw=1, ls='--')
    # Need to add drift as well if present
    ax.set_ylabel('$\Delta$ RV [m/s]')
    if legend:
        ax.legend()
    ax.grid(alpha=0.3)
    if save_fig and fig is not None:
        fig.savefig(f"{star_name}_rv_keplerian_fit.png", dpi=300)
    if show_plot:
        plt.show(block=False)
    return fig, ax


def plot_series(df, quantity, ax=None, fig=None, star_name=None, save_fig=False, show_plot=False, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10, 6)))
    for ins in df['ins_name'].unique():
        ins_mask = df['ins_name'] == ins
        color = ins_colors.get(ins, ins_colors['default'])
        ax.errorbar(df.loc[ins_mask, 'obj_date_bjd'],
                    df.loc[ins_mask, quantity],
                    yerr=df.loc[ins_mask,
                                f"{quantity}_err"] if f"{quantity}_err" in df.columns else None,
                    fmt='o', label=ins, color=color, alpha=0.7,
                    capsize=3, **kwargs)
    ax.set_xlabel('Time [d]')
    ax.set_ylabel(quantity)
    ax.legend()
    ax.grid(alpha=0.3)
    if save_fig and fig is not None:
        fig.savefig(f"{star_name}_{quantity}.png", dpi=300)
    if show_plot:
        plt.show(block=False)
    return fig, ax


def plot_phase_fold(
    rv_model,
    planet_index=0,
    fig=None,
    ax=None,
    star_name=None,
    instruments=None,
    phase_pad=40,
    n_model=1000,
    data_alpha_main=0.7,
    data_alpha_wrap=0.25,
    model_color="black",
    model_lw=1.5,
):
    """
    Plot the phase-folded RV curve of a single Keplerian planet.

    Parameters
    ----------
    rv_model : object
        Fitted RV model containing keplerian components.
    planet_index : int
        Index of planet to phase fold.
    fig, ax : matplotlib Figure/Axes
        Optional external figure and axes.
    star_name : str
        Optional title.
    instruments : array-like
        Instrument labels aligned with rv_model.t.
    color_dic : dict
        Mapping instrument name -> matplotlib color.
    phase_pad : float
        Degrees shown outside [0, 360] for wrap visualization.
    n_model : int
        Number of points for smooth model curve.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    elif fig is None:
        fig = ax.figure
    kep = rv_model.keplerian[str(planet_index)]
    param = kep.get_param()
    # Orbital params
    p = planet_index
    P = rv_model.get_param(f'kep.{p}.P')
    rv_model.set_keplerian_param(str(planet_index), param=['P', 'K', 'la0', 'sqrtecosw', 'sqrtesinw'])
    # Recover e and omega consistently
    if f'kep.{p}.w' in rv_model.fit_param:
        w = rv_model.get_param(f'kep.{p}.w')
    else:
        sqrt_esinw = rv_model.get_param(f'kep.{p}.sqrtesinw')
        sqrt_ecosw = rv_model.get_param(f'kep.{p}.sqrtecosw')
        w = np.arctan2(sqrt_esinw, sqrt_ecosw)

    # Mean longitude reference
    lambda0 = rv_model.get_param(f'kep.{p}.la0')
    rv_model.set_keplerian_param(str(planet_index), param=param)

    M0 = np.degrees(lambda0 - w)

    t = rv_model.t
    kep = rv_model.keplerian[f"{p}"]

    # Residuals cleaned from covariance structure
    res = rv_model.residuals()
    sig = np.sqrt(rv_model.cov.A)
    res -= rv_model.cov.self_conditional(res)

    # Select primary series
    idx = rv_model.series_index[0]
    res = res[idx]
    sig = sig[idx]
    t = t[idx]

    # Mean anomaly in degrees
    M = (M0 + 360.0 / P * t) % 360.0

    # Add planet signal back (so we isolate this planet)
    res_planet = res + kep.rv(t)

    # Smooth model curve
    M_smooth = np.linspace(0, 360, n_model)
    t_model = (M_smooth - M0) * P / 360.0
    model_rv = kep.rv(t_model)

    # Extend for wrap
    M_ext = np.concatenate([M - 360, M, M + 360])
    rv_ext = np.tile(res_planet, 3)
    sig_ext = np.tile(sig, 3)

    ax.plot(
        M_smooth,
        model_rv,
        color=model_color,
        lw=model_lw,
        zorder=3,
    )
    # left and right smooth model
    ax.plot(
        M_smooth-360,
        model_rv,
        color=model_color,
        lw=1,
        alpha=0.2,
        zorder=2,
    )
    ax.plot(
        M_smooth+360,
        model_rv,
        color=model_color,
        lw=1,
        alpha=0.2,
        zorder=2,
    )

    if instruments is None:
        # Plot as single dataset
        main = (M_ext >= 0) & (M_ext <= 360)
        ax.errorbar(
            M_ext[main],
            rv_ext[main],
            sig_ext[main],
            fmt=".",
            color="tab:blue",
            alpha=data_alpha_main,
            markersize=6,
            rasterized=True,
        )
        ax.errorbar(
            M_ext[~main],
            rv_ext[~main],
            sig_ext[~main],
            fmt=".",
            color="tab:blue",
            alpha=data_alpha_wrap,
            markersize=6,
            rasterized=True,
        )
    else:
        instruments_ext = np.tile(instruments[idx], 3)
        unique_inst = np.unique(instruments_ext)

        for inst in unique_inst:
            sel = instruments_ext == inst
            main = sel & (M_ext >= 0) & (M_ext <= 360)

            color = (
                ins_colors[inst]
                if (ins_colors is not None and inst in ins_colors)
                else None
            )

            ax.errorbar(
                M_ext[main],
                rv_ext[main],
                sig_ext[main],
                fmt=".",
                color=color,
                alpha=data_alpha_main,
                markersize=6,
                rasterized=True,
                label=inst
            )

            ax.errorbar(
                M_ext[sel & ~main],
                rv_ext[sel & ~main],
                sig_ext[sel & ~main],
                fmt=".",
                color=color,
                alpha=data_alpha_wrap,
                markersize=6,
                rasterized=True,
            )

    ax.axhline(0, ls="--", lw=1, color="k", alpha=0.3)
    ax.set_xlim(-phase_pad, 360 + phase_pad)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel("Mean anomaly [deg]")
    ax.set_ylabel(r"$\Delta$ RV [m/s]")
    ax.legend()
    if star_name is not None:
        ax.set_title(f"{star_name} - Planet {planet_index+1} - {P:.0f} days")

    fig.tight_layout()

    return fig, ax


def plot_series_periodogram(rv_data, series, adjust_means=True, ax=None, fig=None, fap_level=1e-3, min_freq=None, max_freq=None, samples=10000):
    from .tseries import gls_periodogram
    from .stats import weighted_mean
    rv_data = rv_data.copy()
    if f'{series}_err' in rv_data.columns:
        series_err = rv_data[f'{series}_err']
    else:
        series_err = None
    if (adjust_means):
        for ins in rv_data.ins_name.unique():
            ins_mask = rv_data.ins_name == ins
            mean_rv, _ = weighted_mean(rv_data[series][ins_mask],
                                       series_err[ins_mask])
            rv_data.loc[ins_mask, series] -= mean_rv
            print(f'Adjusted means for ins {ins} by subtracting {mean_rv:.2f}')
    freq, power, best_period, fap, power_threshold = gls_periodogram(
        rv_data.obj_date_bjd, rv_data[series], series_err)
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(1/freq, power)
    ax.set_xscale('log')
    ax.axvline(best_period, ls='--', c='tab:green', alpha=0.5,
               label=f'Highest peak: {best_period:.1f} days, FAP: {fap:.1e}')
    ax.axhline(power_threshold,
               label=f'FAP {fap_level:.1e}', ls='--', c='k', alpha=0.3)
    ax.set_xlabel('Period [d]')
    ax.set_ylabel('Normalised power')
    ax.set_title(series)
    ax.legend()
    return fig, ax
