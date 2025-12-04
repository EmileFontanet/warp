from matplotlib import pyplot as plt

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
    'CORALIE07-NDRS': 'gold',
    'CORALIE14-NDRS': 'darkblue',
    'CORALIE24-NDRS': 'mediumaquamarine',
    'default': 'tab:cyan'
}


def plot_rv(rv_data, ax=None, fig=None, star_name=None, save_fig=False, show_plot=False, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10, 6)))
    for ins in rv_data['ins_name'].unique():
        ins_mask = rv_data['ins_name'] == ins
        color = ins_colors.get(ins, ins_colors['default'])
        ax.errorbar(rv_data.loc[ins_mask, 'obj_date_bjd'],
                    rv_data.loc[ins_mask, 'spectro_ccf_rv'],
                    yerr=rv_data.loc[ins_mask, 'spectro_ccf_rv_err'],
                    fmt='o', label=ins, color=color, alpha=0.7,
                    capsize=3, **kwargs)
    ax.set_xlabel('Time [d]')
    ax.set_ylabel('Radial Velocity [m/s]')
    ax.legend()
    ax.grid(alpha=0.3)
    if save_fig and fig is not None:
        fig.savefig(f"{star_name}_rv.png", dpi=300)
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
