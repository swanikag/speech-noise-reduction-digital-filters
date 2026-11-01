function sos = butter_sos(N, Wn, type)
%BUTTER_SOS  Butterworth filter as second-order sections (gain folded into section 1).
%   type: 'low' | 'high' | 'bandpass' (Wn = [wl wh] normalised to Nyquist).
if strcmp(type, 'bandpass')
    [z, p, k] = butter(N, Wn);          % two-element Wn -> band-pass in MATLAB and Octave
else
    [z, p, k] = butter(N, Wn, type);
end
[sos, g] = zp2sos(z, p, k);
sos(1, 1:3) = sos(1, 1:3) * g;
end
