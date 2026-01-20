from kepmodel import rv
from spleaf import cov, term
import numpy as np

default_inst_jitter = {
    'CORALIE98': 5.0,
    'CORALIE07': 8.0,
    'CORALIE14': 3.0,
    'CORAVEL': 150.0,
    'HARPS03': 0.75,
    'HARPS15': 0.75,
    'HARPN': 0.75,
    'HIRES': 2.5,
    'CHIRON': 8.0,
    'FEROS': 8.0,
    'default': 0.1

}

offsets_relative_to_cor14 = {
    'CORALIE98': -25,
    'CORALIE07': -25,
    'CORALIE24': -10,
    'CORALIE14': 0,
}


def fit_keplerian(time, rv_data, rv_err, instruments, N_pla=3, n_lin=1, stellar_jitter=0, fap_threshold=1e-3, periods_init=[],
                  fix_cor_offsets=False):
    instjit = {}
    for inst in instruments.unique():
        sig = default_inst_jitter.get(
            inst, default_inst_jitter['default'])
        instjit[f'{inst}_jit'] = term.InstrumentJitter(
            instruments == inst, sig)
    instjit['global_jitter'] = term.Jitter(stellar_jitter)
    ref_epoch = np.median(time)
    rv_model = rv.RvModel(
        time - ref_epoch,
        rv_data,
        err=term.Error(rv_err),
        **instjit,

    )
    if (fix_cor_offsets == False or 'CORALIE14' not in instruments.values):
        for inst in instruments.unique():
            rv_model.add_lin(instruments == inst,
                             f'{inst}_offset',
                             value=np.mean(rv_data[instruments == inst]),)
    else:
        rv_model.fixed_offsets = {}
        COR_mask = np.array(['COR' in c for c in instruments])
        # We add a global CORALIE offset
        rv_model.add_lin(COR_mask,
                         f'rv_CORALIE_offset',
                         value=np.mean(rv_data[COR_mask]))
        for inst in instruments.unique():
            # And then add offsets relative to COR14
            if (inst != 'CORALIE14') and (inst in offsets_relative_to_cor14):
                print(
                    f"Fixing offset of {inst} relative to CORALIE14: {offsets_relative_to_cor14[inst]}")
                rv_model.add_lin(instruments == inst,
                                 f'rv_{inst}_offset',
                                 fit=False,
                                 value=offsets_relative_to_cor14[inst])
                rv_model.fixed_offsets[inst] = offsets_relative_to_cor14[inst]
            elif inst != 'CORALIE14':
                rv_model.add_lin(instruments == inst,
                                 f'rv_{inst}_offset',
                                 value=np.mean(rv_data[instruments == inst]))
    rv_model.fit()
    for kpow in range(n_lin):
        rv_model.add_lin(rv_model.t ** (kpow + 1),
                         f"drift_pow{kpow+1}", value=1)
    rv_model.fit()

    Pmax = 1.5*(np.max(time) - np.min(time))
    Pmin = 1.0
    nu0 = 2 * np.pi / Pmax
    nfreq = 50000
    dnu = (2 * np.pi / Pmin - nu0) / (nfreq - 1)
    for p in periods_init:
        rv_model.add_keplerian_from_period(p)
        rv_model.set_keplerian_param(
            f"{rv_model.nkep-1}", param=["P", "la0", "K", "sqrtesinw", "sqrtecosw"]
        )
        rv_model.fit()
    for kpla in range(N_pla-len(periods_init)):
        # Compute periodogram
        nu, power = rv_model.periodogram(nu0, dnu, nfreq)
        P = 2 * np.pi / nu
        # Compute FAP

        kmax = np.argmax(power[P > 2])
        faplvl = rv_model.fap(power[kmax], nu.max())

        if faplvl > fap_threshold:
            break
        # Add new planet
        rv_model.add_keplerian_from_period(P[kmax])
        rv_model.set_keplerian_param(
            f"{rv_model.nkep-1}", param=["P", "la0", "K", "sqrtesinw", "sqrtecosw"]
        )
        # Global fit of the model
        rv_model.fit()
    return rv_model
