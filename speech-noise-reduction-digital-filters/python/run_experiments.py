"""Run all experiments; writes results/*.csv and results/audio/*.wav. Deterministic (seeded)."""
import os, time, json
import numpy as np, pandas as pd
from scipy import signal
from scipy.io import wavfile
from common import *

NSEEDS = 10
out = os.path.join(ROOT, 'results'); os.makedirs(os.path.join(out, 'audio'), exist_ok=True)
clips = load_clips()
filters = design_filters()

# ---------------------------------------------------------------- main trials
rows = []
for ci, (clab, kind, snr_in) in enumerate(CONDITIONS):
    for ki, (cname, s) in enumerate(clips.items()):
        # filtered clean reference per filter (for delay estimation)
        hs_c = {f: causal(sos, s) for f, sos in filters.items()}
        lags = {f: estimate_lag(s, hs_c[f]) for f in filters}
        for seed in range(NSEEDS):
            rng = np.random.default_rng([ci, ki, seed])
            n = scale_noise(s, make_noise(kind, len(s), rng), snr_in)
            x = s + n
            snr0 = snr_db(s, x); mse0 = mse(s, x); mag0 = snr_mag_db(s, x)
            for f, sos in filters.items():
                yc = causal(sos, x)
                yz = zerophase(sos, x)
                ya = align(yc, lags[f])
                oref = {'causal': snr_oref_db(hs_c[f], causal(sos, n)),
                        'zp': snr_oref_db(zerophase(sos, s), zerophase(sos, n))}
                for mode, y, ok in (('causal_raw', yc, 'causal'), ('causal_aligned', ya, 'causal'), ('zero_phase', yz, 'zp')):
                    so = snr_db(s, y)
                    rows.append(dict(condition=clab, clip=cname, seed=seed, filter=f, mode=mode,
                                     snr_in=snr0, snr_out=so, dsnr=so-snr0,
                                     mse_in=mse0, mse=mse(s, y), lag=lags[f],
                                     dsnr_oref=oref[ok]-snr0, dsnr_mag=snr_mag_db(s, y)-mag0))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(out, 'metrics_raw.csv'), index=False)
g = df.groupby(['condition', 'filter', 'mode'], sort=False)
summ = g.agg(snr_in=('snr_in', 'mean'), snr_out=('snr_out', 'mean'), snr_out_std=('snr_out', 'std'),
             dsnr=('dsnr', 'mean'), dsnr_std=('dsnr', 'std'),
             mse_in=('mse_in', 'mean'), mse=('mse', 'mean'), mse_std=('mse', 'std'),
             lag_mean=('lag', 'mean'), dsnr_oref=('dsnr_oref', 'mean'), dsnr_oref_std=('dsnr_oref', 'std'),
             dsnr_mag=('dsnr_mag', 'mean'), dsnr_mag_std=('dsnr_mag', 'std'), n=('dsnr', 'size')).reset_index()
summ.to_csv(os.path.join(out, 'metrics_summary.csv'), index=False)

# ---------------------------------------------------------------- filter response / group delay
f_axis = np.linspace(1, 8000, 16000)
gd_rows = []
notch_c = {'Notch': [50], 'Notch-H': [50, 100, 150]}
for f, sos in filters.items():
    w, h = signal.sosfreqz(sos, worN=f_axis, fs=FS)
    ph = np.unwrap(np.angle(h))
    gd = -np.gradient(ph, 2*np.pi*f_axis)
    mask = (f_axis >= 100) & (f_axis <= 3400)
    for c0 in notch_c.get(f, []):          # exclude +-10 Hz around notch centres (group delay undefined at a zero)
        mask &= np.abs(f_axis - c0) > 10
    for fq in (200, 500, 1000, 2000, 3000):
        gd_rows.append(dict(filter=f, freq=fq, gd_ms=1000*np.interp(fq, f_axis, gd)))
    gd_rows.append(dict(filter=f, freq=-1, gd_ms=1000*np.max(gd[mask])))
pd.DataFrame(gd_rows).to_csv(os.path.join(out, 'group_delay.csv'), index=False)

# ---------------------------------------------------------------- real-time streaming test
rt_rows = []
name = list(clips)[0]; s = clips[name]
rng = np.random.default_rng(2026)
x = s + scale_noise(s, make_noise('white', len(s), rng), 5)
for f, sos in filters.items():
    batch = causal(sos, x)
    for frame in (80, 160, 320, 512):          # 5, 10, 20, 32 ms at 16 kHz
        zi = np.zeros((sos.shape[0], 2)); ys = []; times = []
        for i in range(0, len(x) - frame + 1, frame):
            t0 = time.perf_counter()
            yf, zi = signal.sosfilt(sos, x[i:i+frame], zi=zi)
            times.append(time.perf_counter() - t0); ys.append(yf)
        ys = np.concatenate(ys); L = len(ys)
        times = np.array(times)
        rt_rows.append(dict(filter=f, frame=frame, frame_ms=1000*frame/FS,
                            max_abs_diff=np.max(np.abs(ys - batch[:L])),
                            mean_ms=1000*times.mean(), p99_ms=1000*np.percentile(times, 99),
                            max_ms=1000*times.max(), rtf=times.sum()/(L/FS), nframes=len(times)))
pd.DataFrame(rt_rows).to_csv(os.path.join(out, 'realtime_streaming.csv'), index=False)

# ---------------------------------------------------------------- example audio for listening
rng = np.random.default_rng(2026)
nz = scale_noise(s, make_noise('white', len(s), rng), 5); x = s + nz
def w16(fn, y): wavfile.write(os.path.join(out, 'audio', fn), FS, np.int16(np.clip(y, -1, 1)*32767))
w16(f'{name}_clean.wav', s); w16(f'{name}_noisy_WGN5dB.wav', x)
for f, sos in filters.items():
    w16(f'{name}_WGN5dB_{f}_causal.wav', causal(sos, x)); w16(f'{name}_WGN5dB_{f}_zerophase.wav', zerophase(sos, x))

meta = dict(nseeds=NSEEDS, nclips=len(clips), clips={k: len(v)/FS for k, v in clips.items()}, fs=FS,
            trials_per_cell=NSEEDS*len(clips), numpy=np.__version__)
json.dump(meta, open(os.path.join(out, 'meta.json'), 'w'), indent=1)
print('done', df.shape)
