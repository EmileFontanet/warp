accepted_pipelines = {
    'CORALIE98': '3.3',
    'CORALIE07': '3.4',
    'CORALIE14': '3.8',
    'CORALIE24': '3.8',
    'HARPS03': '3.3.6',
    'HARPS15': '3.3.6',
    'NIRPS': '3.3.12',
    'ESPRESSO19': '3.3.10',
    'ESPRESSO18': '3.3.10',

}
ignored_cols = [
    # 'obj_pos',
    'pub_bibcode',
    'pub_ref',
    'db_spectrum',
    'spectro_analysis',
    'ins_adc',
    'photocenter',
    'spectro_cal_thfile',
    'clim_seeing',
    'spectro_ccf_mask_in_sof',
    'ins_opti',
    'th_ar',
    'ins1_adc',
    'ins2_adc',
    'ins3_adc',
    'ins4_adc'
]
exclude_cols_bin_by_night = [
    'flux',
    'texp',
    'drift',
    'berv',
    'airmass',
    'continuum'
]
