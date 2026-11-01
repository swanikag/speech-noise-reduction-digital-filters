function sos = design_notch_sos(freqs, Q, fs)
%DESIGN_NOTCH_SOS  Cascade of 2nd-order IIR notch filters (same formula as iirnotch).
%   sos = design_notch_sos([50 100 150], 30, 16000) returns one biquad row per notch
%   frequency: [b0 b1 b2 1 a1 a2].
sos = zeros(numel(freqs), 6);
for i = 1:numel(freqs)
    w0   = 2*pi*freqs(i)/fs;
    bw   = w0/Q;
    beta = tan(bw/2);
    g    = 1/(1+beta);
    sos(i,:) = [g*[1, -2*cos(w0), 1], 1, -2*g*cos(w0), 2*g-1];
end
end
