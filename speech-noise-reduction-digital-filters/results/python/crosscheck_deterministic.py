"""RNG-free cross-check (hum with zero phases) mirroring matlab/crosscheck_deterministic.m"""
import os, numpy as np, pandas as pd
from common import *
clips = load_clips(); s = list(clips.values())[0]; F = design_filters()
N = len(s); t = np.arange(N)/FS
n = sum(a*np.sin(2*np.pi*50*(k+1)*t) for k, a in enumerate([1.0, 0.5, 0.25]))
n = scale_noise(s, n, 5); x = s + n; snr0 = snr_db(s, x); mag0 = snr_mag_db(s, x)
rows = []
for f, sos in F.items():
    hs = causal(sos, s); lag = estimate_lag(s, hs)
    yc = causal(sos, x); yz = zerophase(sos, x); ya = align(yc, lag)
    oc = snr_oref_db(hs, causal(sos, n)) - snr0
    oz = snr_oref_db(zerophase(sos, s), zerophase(sos, n)) - snr0
    for mode, y, o in (('causal_raw', yc, oc), ('causal_aligned', ya, oc), ('zero_phase', yz, oz)):
        rows.append(dict(filter=f, mode=mode, snr_out=snr_db(s, y), mse=mse(s, y), dsnr_oref=o, dsnr_mag=snr_mag_db(s, y)-mag0))
pd.DataFrame(rows).to_csv(os.path.join(ROOT, 'results', 'crosscheck_python.csv'), index=False)
m = pd.read_csv(os.path.join(ROOT, 'results', 'crosscheck_matlab.csv'))
p = pd.DataFrame(rows)
d = m.merge(p, on=['filter', 'mode'], suffixes=('_m', '_p'))
for c in ['snr_out', 'mse', 'dsnr_oref', 'dsnr_mag']:
    d['d_'+c] = (d[c+'_m'] - d[c+'_p']).abs()
print(d[['filter', 'mode', 'snr_out_m', 'snr_out_p', 'd_snr_out', 'd_mse', 'd_dsnr_oref', 'd_dsnr_mag']].to_string())
print('MAX abs diff  snr_out %.2e  mse %.2e  oref %.2e  mag %.2e' % tuple(d['d_'+c].max() for c in ['snr_out', 'mse', 'dsnr_oref', 'dsnr_mag']))
