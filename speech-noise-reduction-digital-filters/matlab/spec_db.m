function [D, f, t] = spec_db(x, fs, L, H)
%SPEC_DB  Power spectrogram in dB (periodic Hann, window L, hop H); portable MATLAB/Octave.
x = x(:); nf = floor((numel(x) - L)/H) + 1;
w = 0.5 - 0.5*cos(2*pi*(0:L-1)'/L);
idx = bsxfun(@plus, (1:L)', H*(0:nf-1));
X = fft(bsxfun(@times, x(idx), w), L);
D = 10*log10(abs(X(1:L/2+1, :)).^2 + 1e-14);
f = (0:L/2)' * fs / L; t = ((0:nf-1)*H + L/2) / fs;
end
