function [p, f] = welch_psd(x, nfft, fs)
%WELCH_PSD  Welch PSD with Hann window and 50% overlap (works in MATLAB and Octave).
w = hann(nfft);
if exist('OCTAVE_VERSION', 'builtin')
    [p, f] = pwelch(x, w, 0.5, nfft, fs);      % Octave: overlap as a fraction
else
    [p, f] = pwelch(x, w, nfft/2, nfft, fs);   % MATLAB: overlap in samples
end
end
