%CROSSCHECK_DETERMINISTIC  RNG-free test (hum with fixed zero phases) to compare MATLAB/Octave with Python.
%   Writes ../results/crosscheck_matlab.csv ; compare with python/crosscheck_deterministic.py
clear; if exist('OCTAVE_VERSION', 'builtin'), pkg load signal; end
here = fileparts(mfilename('fullpath')); addpath(here);
fs = 16000; [clips, names] = load_clips(fullfile(here, '..', 'data')); s = clips{1};
F = design_filters(fs); M = metrics();
n = scale_noise(s, make_noise('hum', numel(s), fs, [0 0 0]), 5); x = s + n; snr0 = M.snr_db(s, x);
fid = fopen(fullfile(here, '..', 'results', 'crosscheck_matlab.csv'), 'w');
fprintf(fid, 'filter,mode,snr_out,mse,dsnr_oref,dsnr_mag\n');
for fi = 1:numel(F)
    sos = F(fi).sos; hs = sos_filter(sos, s); lag = M.lag(s, hs, 160);
    yc = sos_filter(sos, x); yz = zero_phase(sos, x);
    oc = M.snr_oref(hs, sos_filter(sos, n)) - snr0;
    oz = M.snr_oref(zero_phase(sos, s), zero_phase(sos, n)) - snr0;
    ya = M.align(yc, lag);
    fprintf(fid, '%s,causal_raw,%.6f,%.8e,%.6f,%.6f\n', F(fi).name, M.snr_db(s,yc), M.mse(s,yc), oc, M.snr_mag(s,yc)-M.snr_mag(s,x));
    fprintf(fid, '%s,causal_aligned,%.6f,%.8e,%.6f,%.6f\n', F(fi).name, M.snr_db(s,ya), M.mse(s,ya), oc, M.snr_mag(s,ya)-M.snr_mag(s,x));
    fprintf(fid, '%s,zero_phase,%.6f,%.8e,%.6f,%.6f\n', F(fi).name, M.snr_db(s,yz), M.mse(s,yz), oz, M.snr_mag(s,yz)-M.snr_mag(s,x));
end
fclose(fid); disp('wrote results/crosscheck_matlab.csv');
