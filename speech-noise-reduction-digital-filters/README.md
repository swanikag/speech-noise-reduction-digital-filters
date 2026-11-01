# Real-Time Speech Noise Reduction Using Digital Filters

Code, data, results and IEEE paper source comparing **Butterworth LPF, HPF, BPF and IIR notch filters** for speech
noise reduction, evaluated with SNR, SNR improvement (dSNR) and MSE, in both **causal (real-time, frame-by-frame)**
and **offline zero-phase** form.

## Key findings (measured, 6 clips x 7 noise conditions x 10 seeds)
| Condition (5 dB input) | Best causal filter | dSNR (dB) | Zero-phase |
|---|---|---|---|
| White noise | LPF | +2.89 | LPF +3.28, BPF +3.15 |
| 50 Hz hum + harmonics | Notch bank (50/100/150 Hz) | +11.49 | +16.70 |
| High-frequency hiss | LPF | +11.46 | - |

Causal HPF/BPF lose waveform SNR (group delay ~6 ms near 100 Hz) although they remove noise power
(output-referred dSNR up to 13.6 dB for rumble). Streaming output equals batch causal output exactly.

## Repository layout
```
matlab/    MATLAB/Octave implementation (run_experiments.m, streaming_demo.m, plot_figures.m, ...) + matlab_code.pdf
python/    Pipeline that produced the reported numbers and figures
data/      6 speech clips (see data/README.md for licence notes)
results/   per-trial CSVs, summaries, group delay, streaming timing, example WAVs
figures/   all figures (PDF + 600-dpi PNG) and TikZ sources
paper/     IEEEtran LaTeX source (main.tex, references.bib), compiled main.pdf
```

## Reproduce
```bash
pip install -r requirements.txt
cd python
python run_experiments.py      # metrics CSVs (seeded, deterministic)
python make_figures.py         # figures -> ../figures
python make_tables.py          # LaTeX tables -> ../results/tab_*.tex
python crosscheck_deterministic.py   # compares with Octave/MATLAB output
cp ../results/tab_*.tex ../paper/ && cp ../figures/*.pdf ../paper/figures/
cd ../paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
```
MATLAB / GNU Octave (with `signal` package): `cd matlab; run_experiments` (and `plot_figures`).

## Validation and honest notes
- Reported results come from the Python code. The MATLAB code was executed in **GNU Octave 8.4.0** (not MATLAB itself) and
  matches Python to ~5e-7 dB on a deterministic test and within 0.12 dB (mean 0.004 dB) on randomised runs.
- `matlab/plot_figures.m` could not be rendered in the authors' environment and is untested.
- Noise is synthetic; only 6 clips from 2 speakers; no PESQ/STOI or listening tests. Timings were measured in a Python sandbox.

## Paper
Overleaf: upload this repo's `paper/` folder (or the separate Overleaf ZIP). Replace the placeholder author/affiliation in `paper/main.tex`.

## License
Code: MIT. Speech data: see `data/README.md`.
