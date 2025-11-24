from astropy.timeseries import LombScargle
import numpy as np


def gls_periodogram(t, y, yerr=None, min_freq=None, max_freq=None, samples=10000):

    ls = LombScargle(t, y, yerr, center_data=True)
    if min_freq is None:
        min_freq = 1 / (max(t) - min(t))
    if max_freq is None:
        dt = np.median(np.diff(np.sort(t)))
        max_freq = 0.5 / dt
        max_freq = max(max_freq, 1.0)
    freq = np.linspace(min_freq, max_freq, samples)
    power = ls.power(freq)

    best_freq = freq[np.argmax(power)]
    best_period = 1 / best_freq
    fap = ls.false_alarm_probability(power.max())
    return freq, power, best_period, fap
