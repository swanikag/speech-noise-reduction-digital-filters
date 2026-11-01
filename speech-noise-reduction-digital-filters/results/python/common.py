"""Shared helpers: filter design, noise generation, metrics (mirrors matlab/*.m)."""
import numpy as np
from scipy import signal
from scipy.io import wavfile
import glob, os

FS = 16000
TARGET_RMS = 0.1          # every clean clip is scaled to RMS = 0.1 (-20 dBFS)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_clips(folder=None):
    folder = folder or os.path.join(ROOT, 'data')
    clips = {}
    for f in sorted(glob.glob(os.path.join(folder, '*.wav'))):
        fs, x = wavfile.read(f)
        assert fs == FS, f"{f}: expected {FS} Hz"
        x = x.astype(np.float64) / 32768.0
        x = x - np.mean(x)
        x = x * (TARGET_RMS / np.sqrt(np.mean(x**2)))
        clips[os.path.basename(f)[:-4]] = x
    return clips


def notch_sos(freqs, Q=30.0, fs=FS):
    """Cascade of second-order IIR notches (same formula as MATLAB iirnotch)."""
    rows = []
    for f0 in freqs:
        w0 = 2*np.pi*f0/fs
        bw = w0/Q
        beta = np.tan(bw/2)
        g = 1.0/(1.0+beta)
        b = g*np.array([1.0, -2*np.cos(w0), 1.0])
        a = np.array([1.0, -2*g*np.cos(w0), 2*g-1.0])
        rows.append(np.concatenate([b, a]))
    return np.array(rows)


def design_filters(fs=FS):
    return {
        'LPF':     signal.butter(4, 3400, 'low', fs=fs, output='sos'),
        'HPF':     signal.butter(4, 100, 'high', fs=fs, output='sos'),
        'BPF':     signal.butter(4, [100, 3400], 'bandpass', fs=fs, output='sos'),
        'Notch':   notch_sos([50]),
        'Notch-H': notch_sos([50, 100, 150]),
    }


def causal(sos, x):
    return signal.sosfilt(sos, x)


def zerophase(sos, x):
    """Offline forward-backward filtering (non-causal), no edge padding (same as matlab/zerophase.m)."""
    y = signal.sosfilt(sos, x)
    y = signal.sosfilt(sos, y[::-1])
    return y[::-1]


def scale_noise(s, n, snr_db):
    return n*np.sqrt(np.sum(s**2)/(np.sum(n**2)*10**(snr_db/10)))


def make_noise(kind, N, rng, fs=FS):
    t = np.arange(N)/fs
    if kind == 'white':
        return rng.standard_normal(N)
    if kind == 'hum':
        ph = rng.uniform(0, 2*np.pi, 3)
        return sum(a*np.sin(2*np.pi*50*(k+1)*t+ph[k]) for k, a in enumerate([1.0, 0.5, 0.25]))
    if kind == 'rumble':
        sos = signal.butter(4, [20, 80], 'bandpass', fs=fs, output='sos')
        return signal.sosfilt(sos, rng.standard_normal(N))
    if kind == 'hiss':
        sos = signal.butter(4, [5000, 7900], 'bandpass', fs=fs, output='sos')
        return signal.sosfilt(sos, rng.standard_normal(N))
    if kind == 'hum+white':
        h = make_noise('hum', N, rng); w = make_noise('white', N, rng)
        return h/np.sqrt(np.mean(h**2)) + w/np.sqrt(np.mean(w**2))
    raise ValueError(kind)


CONDITIONS = [  # (label, noise kind, input SNR dB)
    ('WGN 0 dB',    'white',     0),
    ('WGN 5 dB',    'white',     5),
    ('WGN 10 dB',   'white',    10),
    ('Hum 5 dB',    'hum',       5),
    ('Rumble 5 dB', 'rumble',    5),
    ('Hiss 5 dB',   'hiss',      5),
    ('Hum+WGN 5 dB','hum+white', 5),
]


def snr_db(s, y):
    return 10*np.log10(np.sum(s**2)/np.sum((s-y)**2))


def mse(s, y):
    return np.mean((s-y)**2)


def estimate_lag(s, hs, maxlag=160):
    """Integer lag (samples) maximizing cross-correlation between clean s and filtered clean hs."""
    c = signal.correlate(hs, s, mode='full', method='fft')
    lags = signal.correlation_lags(len(hs), len(s), mode='full')
    m = (lags >= 0) & (lags <= maxlag)
    return int(lags[m][np.argmax(c[m])])


def align(y, lag):
    """Advance y by lag samples (compensate delay); zero-pad tail to keep length."""
    if lag == 0:
        return y
    return np.concatenate([y[lag:], np.zeros(lag)])


def stft_mag(x, fs=FS):
    _, _, Z = signal.stft(x, fs=fs, window='hann', nperseg=512, noverlap=256, boundary=None, padded=False)
    return np.abs(Z)


def snr_mag_db(s, y):
    """Phase-insensitive SNR on STFT magnitudes (512-pt Hann, 50% overlap)."""
    S, Y = stft_mag(s), stft_mag(y)
    return 10*np.log10(np.sum(S**2)/np.sum((S-Y)**2))


def snr_oref_db(hs, hn):
    """Output-referred SNR: power of filtered clean speech over power of filtered noise (delay invariant)."""
    return 10*np.log10(np.sum(hs**2)/np.sum(hn**2))
