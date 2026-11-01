"""Generate all result figures (PDF vector + 600-dpi PNG) into ../figures from ../results CSVs."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal
from common import *

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['STIXGeneral', 'DejaVu Serif'], 'mathtext.fontset': 'stix',
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7,
    'ytick.labelsize': 7, 'axes.linewidth': 0.6, 'lines.linewidth': 1.0, 'grid.linewidth': 0.4,
    'grid.alpha': 0.5, 'savefig.dpi': 600, 'figure.dpi': 150, 'pdf.fonttype': 42})
W1, W2 = 3.5, 7.16
COL = {'LPF': '#0072B2', 'HPF': '#D55E00', 'BPF': '#009E73', 'Notch': '#CC79A7', 'Notch-H': '#E69F00'}
FIG = os.path.join(ROOT, 'figures'); os.makedirs(FIG, exist_ok=True)
F = design_filters()
summ = pd.read_csv(os.path.join(ROOT, 'results', 'metrics_summary.csv'))
cond_order = [c[0] for c in CONDITIONS]


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + '.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(FIG, name + '.png'), bbox_inches='tight')
    plt.close(fig)


def panel(ax, txt):
    ax.set_title(txt, loc='left', fontsize=8, fontweight='bold', pad=2)


# ------------------------------------------------------------ Fig: filter responses
fq = np.logspace(np.log10(20), np.log10(8000), 4000)
fig, ax = plt.subplots(1, 3, figsize=(W2, 2.3))
for k, sos in F.items():
    _, h = signal.sosfreqz(sos, worN=fq, fs=FS)
    ax[0].semilogx(fq, 20*np.log10(np.maximum(abs(h), 1e-6)), color=COL[k], label=k,
                   ls={'LPF': '-', 'HPF': '--', 'BPF': '-.', 'Notch': '-', 'Notch-H': ':'}[k],
                   lw={'BPF': 1.6, 'Notch-H': 1.4}.get(k, 1.0))
ax[0].set_ylim(-60, 5); ax[0].set_xlim(20, 8000); ax[0].set_xlabel('Frequency (Hz)'); ax[0].set_ylabel('Magnitude (dB)')
ax[0].grid(True, which='both'); ax[0].legend(loc='center', bbox_to_anchor=(0.55, 0.42), ncol=1, frameon=True, framealpha=0.9, edgecolor='none')
panel(ax[0], '(a)')
for k, ls in (('LPF', '-'), ('HPF', '-')):
    _, h = signal.sosfreqz(F[k], worN=fq, fs=FS)
    ax[1].semilogx(fq, 20*np.log10(np.maximum(abs(h), 1e-6)), ls, color=COL[k], label=f'{k} causal')
    ax[1].semilogx(fq, 40*np.log10(np.maximum(abs(h), 1e-6)), '--', color=COL[k], label=f'{k} zero-phase')
ax[1].axhline(-3.01, color='k', lw=0.4, ls=':'); ax[1].axhline(-6.02, color='k', lw=0.4, ls=':')
ax[1].set_ylim(-60, 5); ax[1].set_xlim(20, 8000); ax[1].set_xlabel('Frequency (Hz)'); ax[1].set_ylabel('Magnitude (dB)')
ax[1].grid(True, which='both'); ax[1].legend(loc='center', bbox_to_anchor=(0.55, 0.42), ncol=1, frameon=True, framealpha=0.9, edgecolor='none'); panel(ax[1], '(b)')
fl = np.linspace(20, 4000, 40000)
notch_c = {'Notch': [50], 'Notch-H': [50, 100, 150]}
for k, sos in F.items():
    _, h = signal.sosfreqz(sos, worN=fl, fs=FS)
    gd = -np.gradient(np.unwrap(np.angle(h)), 2*np.pi*fl)*1000
    for c0 in notch_c.get(k, []):
        gd[np.abs(fl - c0) < 6] = np.nan
    ax[2].semilogx(fl, gd, color=COL[k], label=k)
ax[2].set_ylim(0, 8); ax[2].set_xlim(20, 4000); ax[2].set_xlabel('Frequency (Hz)'); ax[2].set_ylabel('Group delay (ms)')
ax[2].grid(True, which='both'); panel(ax[2], '(c)')
fig.tight_layout(); save(fig, 'fig_filter_response')

# ------------------------------------------------------------ example trial (WGN 5 dB, first clip)
clips = load_clips(); name = list(clips)[0]; s = clips[name]
rng = np.random.default_rng(2026)
n = scale_noise(s, make_noise('white', len(s), rng), 5); x = s + n
lags = {k: estimate_lag(s, causal(F[k], s)) for k in F}
y = {k: align(causal(F[k], x), lags[k]) for k in F}

# Fig: waveforms
t = np.arange(len(s))/FS
fig, axs = plt.subplots(4, 1, figsize=(W1, 4.2), sharex=True, sharey=True)
for a, (sig, lab) in zip(axs, ((s, 'Clean'), (x, 'Noisy (WGN, 5 dB)'), (y['LPF'], 'LPF (causal)'), (y['BPF'], 'BPF (causal)'))):
    a.plot(t, sig, lw=0.4, color='k'); a.set_ylabel('Amplitude'); a.set_ylim(-0.6, 0.6)
    a.text(0.99, 0.93, lab, transform=a.transAxes, ha='right', va='top', fontsize=7)
    a.grid(True)
axs[-1].set_xlabel('Time (s)'); axs[-1].set_xlim(0, 3.0)
fig.tight_layout(); save(fig, 'fig_waveforms')

# Fig: spectra (Welch PSD)
def psd(sig, nper=1024):
    f, p = signal.welch(sig, FS, window='hann', nperseg=nper, noverlap=nper//2)
    return f, 10*np.log10(p + 1e-20)
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6))
for sig, lab, c, ls in ((s, 'Clean', 'k', '-'), (x, 'Noisy', '0.55', '-'), (y['LPF'], 'LPF', COL['LPF'], '-'), (y['BPF'], 'BPF', COL['BPF'], '--')):
    f, p = psd(sig); ax[0].plot(f, p, color=c, ls=ls, label=lab, lw=0.9)
ax[0].set_xlabel('Frequency (Hz)'); ax[0].set_ylabel('PSD (dB/Hz)'); ax[0].set_xlim(0, 8000); ax[0].grid(True)
ax[0].legend(frameon=False, loc='lower left'); panel(ax[0], '(a)')
rng2 = np.random.default_rng(7); nh = scale_noise(s, make_noise('hum', len(s), rng2), 5); xh = s + nh
for sig, lab, c, ls in ((xh, 'Noisy (hum, 5 dB)', '0.55', '-'),
                        (causal(F['Notch'], xh), 'Notch (50 Hz)', COL['Notch'], '-'),
                        (causal(F['Notch-H'], xh), 'Notch-H (50/100/150 Hz)', COL['Notch-H'], '-'),
                        (causal(F['HPF'], xh), 'HPF (100 Hz)', COL['HPF'], '--')):
    f, p = psd(sig, 4096); m = f <= 400; ax[1].plot(f[m], p[m], color=c, ls=ls, label=lab, lw=0.9)
for c0 in (50, 100, 150):
    ax[1].axvline(c0, color='k', lw=0.3, ls=':')
ax[1].set_xlabel('Frequency (Hz)'); ax[1].set_ylabel('PSD (dB/Hz)'); ax[1].set_xlim(0, 400); ax[1].grid(True)
ax[1].legend(frameon=False, loc='lower right', fontsize=6); panel(ax[1], '(b)')
fig.tight_layout(); save(fig, 'fig_spectra')

# Fig: spectrograms
fig, axs = plt.subplots(1, 3, figsize=(W2, 2.2), sharey=True)
vmax = None
for a, (sig, lab) in zip(axs, ((s, 'Clean'), (x, 'Noisy (WGN, 5 dB)'), (y['BPF'], 'BPF output (causal)'))):
    f, tt, S = signal.spectrogram(sig, FS, window='hann', nperseg=512, noverlap=384)
    D = 10*np.log10(S + 1e-14); vmax = D.max() if vmax is None else vmax
    im = a.pcolormesh(tt, f/1000, D, vmin=vmax-70, vmax=vmax, shading='auto', cmap='magma', rasterized=True)
    a.set_title(lab); a.set_xlabel('Time (s)')
axs[0].set_ylabel('Frequency (kHz)')
cb = fig.colorbar(im, ax=axs, pad=0.01, fraction=0.025); cb.set_label('dB')
save(fig, 'fig_spectrograms')

# ------------------------------------------------------------ Fig: dSNR bars
def bars(ax, key, mode, std_key=None, ylabel='', show_leg=False):
    d = summ[summ['mode'] == mode]; wd = 0.16
    for i, k in enumerate(F):
        v = d[d['filter'] == k].set_index('condition').loc[cond_order]
        e = v[std_key] if std_key else None
        ax.bar(np.arange(len(cond_order)) + (i-2)*wd, v[key], wd, yerr=e, color=COL[k], label=k,
               error_kw=dict(lw=0.5, capsize=0.8), linewidth=0)
    ax.set_xticks(np.arange(len(cond_order))); ax.set_xticklabels(cond_order, rotation=30, ha='right')
    ax.axhline(0, color='k', lw=0.5); ax.set_ylabel(ylabel); ax.grid(True, axis='y')
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.7))
bars(ax[0], 'dsnr', 'causal_aligned', 'dsnr_std', r'$\Delta$SNR, waveform (dB)'); panel(ax[0], '(a)')
bars(ax[1], 'dsnr_oref', 'causal_aligned', 'dsnr_oref_std', r'$\Delta$SNR, output-referred (dB)'); panel(ax[1], '(b)')
ax[1].legend(frameon=False, ncol=2, loc='upper left', fontsize=6)
fig.tight_layout(); save(fig, 'fig_dsnr_causal')

# Fig: causal raw vs aligned vs zero-phase
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6)); wd = 0.26
for a, k, lab in ((ax[0], 'BPF', '(a) BPF'), (ax[1], 'LPF', '(b) LPF')):
    for i, (mode, c, ml) in enumerate((('causal_raw', '#999999', 'causal, no delay comp.'),
                                       ('causal_aligned', COL[k], 'causal, delay-compensated'),
                                       ('zero_phase', '#222222', 'zero-phase (offline)'))):
        v = summ[(summ['mode'] == mode) & (summ['filter'] == k)].set_index('condition').loc[cond_order]
        a.bar(np.arange(len(cond_order)) + (i-1)*wd, v['dsnr'], wd, yerr=v['dsnr_std'], color=c, label=ml,
              error_kw=dict(lw=0.5, capsize=0.8), linewidth=0)
    a.set_xticks(np.arange(len(cond_order))); a.set_xticklabels(cond_order, rotation=30, ha='right')
    a.axhline(0, color='k', lw=0.5); a.set_ylabel(r'$\Delta$SNR, waveform (dB)'); a.grid(True, axis='y')
    a.text(0.02, 0.97, lab, transform=a.transAxes, va='top', fontsize=8, fontweight='bold')
ax[0].legend(frameon=False, fontsize=6, loc='upper right')
fig.tight_layout(); save(fig, 'fig_causal_vs_zerophase')

# Fig: MSE (log) causal-aligned vs zero-phase, WGN 5 dB + Hum 5 dB
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.4)); wd = 0.35
for a, cond in ((ax[0], 'WGN 5 dB'), (ax[1], 'Hum 5 dB')):
    base = summ[(summ['condition'] == cond)]['mse_in'].iloc[0]
    for i, (mode, ml, c) in enumerate((('causal_aligned', 'causal (delay-comp.)', '#4C72B0'), ('zero_phase', 'zero-phase', '#222222'))):
        v = summ[(summ['mode'] == mode) & (summ['condition'] == cond)].set_index('filter').loc[list(F)]
        a.bar(np.arange(5) + (i-0.5)*wd, v['mse'], wd, yerr=v['mse_std'], color=c, label=ml, error_kw=dict(lw=0.5, capsize=0.8), linewidth=0)
    a.axhline(base, color='r', lw=0.8, ls='--', label='unfiltered input')
    a.set_yscale('log'); a.set_xticks(range(5)); a.set_xticklabels(list(F)); a.set_ylabel('MSE'); a.grid(True, axis='y', which='both')
    a.set_title(cond)
ax[0].legend(frameon=False, fontsize=6, loc='upper left')
fig.tight_layout(); save(fig, 'fig_mse')

# Fig: real-time budget
rt = pd.read_csv(os.path.join(ROOT, 'results', 'realtime_streaming.csv'))
gd = pd.read_csv(os.path.join(ROOT, 'results', 'group_delay.csv'))
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.4))
for k in F:
    v = rt[rt['filter'] == k]
    ax[0].plot(v['frame_ms'], v['max_ms'], 'o-', ms=3, color=COL[k], label=k)
xx = np.array([4, 40]); ax[0].plot(xx, xx, 'k--', lw=0.7, label='deadline (= frame length)')
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('Frame length (ms)'); ax[0].set_ylabel('Worst-case time per frame (ms)')
ax[0].grid(True, which='both'); ax[0].legend(frameon=False, fontsize=6, loc='lower right'); panel(ax[0], '(a)')
frame_ms = 10.0; names = list(F)
g = gd[gd['freq'] == -1].set_index('filter').loc[names]['gd_ms']
ax[1].bar(range(5), [frame_ms]*5, color='#BBBBBB', label='buffering (10 ms frame)', linewidth=0)
ax[1].bar(range(5), g, bottom=[frame_ms]*5, color=[COL[k] for k in names], label='max group delay, 100-3400 Hz', linewidth=0)
ax[1].set_xticks(range(5)); ax[1].set_xticklabels(names); ax[1].set_ylabel('Algorithmic latency (ms)')
ax[1].grid(True, axis='y'); ax[1].legend(frameon=False, fontsize=6, loc='upper right'); panel(ax[1], '(b)')
ax[1].set_ylim(0, 22)
fig.tight_layout(); save(fig, 'fig_realtime')
print('figures written:', sorted(os.listdir(FIG)))
