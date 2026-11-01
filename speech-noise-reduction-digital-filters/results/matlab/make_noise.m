function n = make_noise(kind, N, fs, phases)
%MAKE_NOISE  Synthetic noise types used in the experiments (unscaled).
%   kind: 'white' | 'hum' | 'rumble' | 'hiss' | 'hum+white'
%   phases (optional, 1x3): fixed hum phases for deterministic tests.
t = (0:N-1)'/fs;
switch kind
    case 'white'
        n = randn(N, 1);
    case 'hum'
        if nargin < 4, phases = 2*pi*rand(1, 3); end
        n = zeros(N, 1); amps = [1.0 0.5 0.25];
        for k = 1:3
            n = n + amps(k)*sin(2*pi*50*k*t + phases(k));
        end
    case 'rumble'
        n = sos_filter(butter_sos(4, [20 80]/(fs/2), 'bandpass'), randn(N, 1));
    case 'hiss'
        n = sos_filter(butter_sos(4, [5000 7900]/(fs/2), 'bandpass'), randn(N, 1));
    case 'hum+white'
        h = make_noise('hum', N, fs); w = make_noise('white', N, fs);
        n = h/sqrt(mean(h.^2)) + w/sqrt(mean(w.^2));
    otherwise
        error('unknown noise kind: %s', kind);
end
end
