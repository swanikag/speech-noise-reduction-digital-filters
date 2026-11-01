"""Create LaTeX tables (results/tables.tex -> \\input in main.tex) strictly from results/*.csv"""
import os, numpy as np, pandas as pd
from common import *
S = pd.read_csv(os.path.join(ROOT, 'results', 'metrics_summary.csv'))
RT = pd.read_csv(os.path.join(ROOT, 'results', 'realtime_streaming.csv'))
GD = pd.read_csv(os.path.join(ROOT, 'results', 'group_delay.csv'))
FILT = ['LPF', 'HPF', 'BPF', 'Notch', 'Notch-H']; CONDS = [c[0] for c in CONDITIONS]
def fmt(x, dec):
    s = f"{x:.{dec}f}"
    if float(s) == 0: s = s.lstrip('-')
    return s.replace('-', '$-$')

def esc(c): return c.replace('+', '+').replace('%', '\\%')

def grid(mode, mean, std, best, scale=1.0, dec=2, extra=None, cap='', lab='', note=''):
    d = S[S['mode'] == mode]
    head = ' & '.join(['Condition'] + (['Unfilt.'] if extra else []) + FILT)
    lines = [r'\begin{table*}[!t]', r'\caption{' + cap + r'}', r'\label{' + lab + r'}', r'\centering\footnotesize',
             r'\begin{tabular}{l' + 'c'*(len(FILT) + (1 if extra else 0)) + '}', r'\hline', head + r' \\ \hline']
    for c in CONDS:
        v = d[d['condition'] == c].set_index('filter').loc[FILT]
        vals = v[mean].values*scale; sds = v[std].values*scale if std else None
        b = np.argmax(vals) if best == 'max' else np.argmin(vals)
        cells = [esc(c)]
        if extra: cells.append(fmt(d[d['condition']==c][extra].iloc[0]*scale, dec))
        for i, f in enumerate(FILT):
            t = fmt(vals[i], dec) + (f" $\\pm$ {sds[i]:.{dec}f}" if std else '')
            cells.append(r'\textbf{' + t + '}' if i == b else t)
        lines.append(' & '.join(cells) + r' \\')
    lines += [r'\hline', r'\end{tabular}', (r'\par\smallskip\raggedright\footnotesize ' + note) if note else '', r'\end{table*}']
    return '\n'.join(lines)

out = []
out.append(grid('causal_aligned', 'dsnr', 'dsnr_std', 'max', cap=r'Waveform $\Delta$SNR (dB), real-time causal filtering with integer delay compensation (mean $\pm$ std over 60 trials; best per row in bold).', lab='tab:dsnr_causal'))
out.append(grid('causal_aligned', 'mse', 'mse_std', 'min', scale=1e3, dec=2, extra='mse_in', cap=r'MSE against the clean reference ($\times10^{-3}$, clips normalised to RMS 0.1), causal filtering with delay compensation (mean $\pm$ std; lowest per row in bold).', lab='tab:mse_causal'))
out.append(grid('zero_phase', 'dsnr', 'dsnr_std', 'max', cap=r'Waveform $\Delta$SNR (dB), offline zero-phase (forward--backward) filtering (mean $\pm$ std over 60 trials; best per row in bold).', lab='tab:dsnr_zp'))
out.append(grid('zero_phase', 'mse', 'mse_std', 'min', scale=1e3, dec=2, extra='mse_in', cap=r'MSE ($\times10^{-3}$), offline zero-phase filtering (mean $\pm$ std; lowest per row in bold).', lab='tab:mse_zp'))
out.append(grid('causal_aligned', 'dsnr_oref', 'dsnr_oref_std', 'max', cap=r'Output-referred $\Delta$SNR (dB), $10\log_{10}(\|Hs\|^2/\|Hn\|^2)-\mathrm{SNR_{in}}$, causal filtering. This metric is independent of filter delay and phase and isolates noise-power removal.', lab='tab:dsnr_oref'))
out.append(grid('causal_aligned', 'dsnr_mag', 'dsnr_mag_std', 'max', cap=r'STFT-magnitude $\Delta$SNR (dB), causal filtering (phase-insensitive; 512-point Hann window, 50\% overlap).', lab='tab:dsnr_mag'))

# real-time table
g = GD[GD['freq'] == -1].set_index('filter')['gd_ms']
g200 = GD[GD['freq'] == 200].set_index('filter')['gd_ms']; g1k = GD[GD['freq'] == 1000].set_index('filter')['gd_ms']
lagm = S[(S['mode'] == 'causal_aligned') & (S['condition'] == 'WGN 5 dB')].set_index('filter')['lag_mean']
r10 = RT[RT['frame'] == 160].set_index('filter')
rows = []
order = {'LPF': 4, 'HPF': 4, 'BPF': 8, 'Notch': 2, 'Notch-H': 6}
for f in FILT:
    rows.append(' & '.join([f, str(order[f]), f"{g200[f]:.2f}", f"{g1k[f]:.2f}", f"{g[f]:.2f}", f"{lagm[f]:.1f} ({lagm[f]/16:.2f})",
                            f"{r10.loc[f,'mean_ms']:.3f}", f"{r10.loc[f,'max_ms']:.3f}", f"{r10.loc[f,'rtf']:.4f}", f"{RT[RT['filter']==f]['max_abs_diff'].max():.1e}"]) + r' \\')
out.append('\n'.join([r'\begin{table*}[!t]', r'\caption{Real-time characteristics of the causal filters at $f_s=16$ kHz (frame-by-frame streaming, 10~ms frames). Group delay (GD) from the analytic phase response; ``max GD'' is the maximum over 100--3400~Hz (excluding $\pm10$~Hz around notch centres, where GD is undefined). Timing was measured on the authors\' sandbox CPU with SciPy and will differ on other hardware.}', r'\label{tab:rt}', r'\centering\footnotesize',
 r'\begin{tabular}{lccccccccc}', r'\hline', r'Filter & Order & GD@200\,Hz & GD@1\,kHz & max GD & est. lag & mean time & max time & RTF & max$|$stream$-$batch$|$ \\',
 r' & & (ms) & (ms) & (ms) & samples (ms) & (ms/frame) & (ms/frame) & & \\ \hline'] + rows + [r'\hline', r'\end{tabular}', r'\end{table*}']))
open(os.path.join(ROOT, 'results', 'tables.tex'), 'w').write('\n\n'.join(out))
print('tables.tex written', sum(len(o) for o in out))

# split into one file per table so main.tex can \input them (copy results/tab_*.tex into paper/)
import re
for p in re.findall(r'\\begin\{table\*\}.*?\\end\{table\*\}', '\n\n'.join(out), re.S):
    lab = re.search(r'\\label\{(.*?)\}', p).group(1)
    open(os.path.join(ROOT, 'results', lab.replace(':', '_') + '.tex'), 'w').write(p)
print('per-table files written to results/tab_*.tex')
