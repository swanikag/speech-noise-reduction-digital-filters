function M = metrics()
%METRICS  Function handles for the evaluation metrics used in the paper.
M.snr_db   = @(s, y) 10*log10(sum(s.^2) / sum((s - y).^2));      % waveform SNR vs clean reference
M.mse      = @(s, y) mean((s - y).^2);                           % mean squared error vs clean
M.snr_oref = @(hs, hn) 10*log10(sum(hs.^2) / sum(hn.^2));        % output-referred (delay invariant)
M.snr_mag  = @snr_mag_db;                                        % STFT-magnitude SNR (phase insensitive)
M.lag      = @estimate_lag;
M.align    = @(y, lag) [y(lag+1:end); zeros(lag, 1)];
end

function v = snr_mag_db(s, y)
S = stft_mag(s); Y = stft_mag(y);
v = 10*log10(sum(S(:).^2) / sum((S(:) - Y(:)).^2));
end

function A = stft_mag(x)
% 512-point periodic Hann, hop 256, no padding, one-sided (257 bins)
x = x(:); L = 512; H = 256; nf = floor((numel(x) - L)/H) + 1;
w = 0.5 - 0.5*cos(2*pi*(0:L-1)'/L);
idx = bsxfun(@plus, (1:L)', H*(0:nf-1));
X = fft(bsxfun(@times, x(idx), w), L);
A = abs(X(1:L/2+1, :));
end

function lag = estimate_lag(s, hs, maxlag)
% integer lag (samples, 0..maxlag) maximising cross-correlation of filtered clean hs with clean s
if nargin < 3, maxlag = 160; end
s = s(:); hs = hs(:); c = zeros(maxlag+1, 1);
for k = 0:maxlag
    c(k+1) = sum(hs(1+k:end) .* s(1:end-k));
end
[~, i] = max(c); lag = i - 1;
end
