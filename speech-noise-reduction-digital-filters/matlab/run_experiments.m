%RUN_EXPERIMENTS  Main experiment: 5 Butterworth/notch filters x 7 noise conditions x clips x seeds.
%   Evaluates causal (real-time) and zero-phase (offline) filtering with waveform SNR, delta-SNR and
%   MSE (+ output-referred and STFT-magnitude delta-SNR). Writes ../results/matlab_metrics_summary.csv
%   Runs unchanged in MATLAB (Signal Processing Toolbox) and GNU Octave (signal package).
clear; clc;
if exist('OCTAVE_VERSION', 'builtin'), pkg load signal; end
here = fileparts(mfilename('fullpath'));
addpath(here);
fs = 16000; NSEEDS = 10;
[clips, names] = load_clips(fullfile(here, '..', 'data'));
F = design_filters(fs); M = metrics();
cond = {'WGN 0 dB','white',0; 'WGN 5 dB','white',5; 'WGN 10 dB','white',10; 'Hum 5 dB','hum',5; ...
        'Rumble 5 dB','rumble',5; 'Hiss 5 dB','hiss',5; 'Hum+WGN 5 dB','hum+white',5};
modes = {'causal_raw','causal_aligned','zero_phase'};
nC = size(cond,1); nF = numel(F); nM = numel(modes); nT = numel(clips)*NSEEDS;
R.snr_out = zeros(nC,nF,nM,nT); R.dsnr = R.snr_out; R.mse = R.snr_out; R.oref = R.snr_out; R.mag = R.snr_out;
for ci = 1:nC
    t = 0;
    for ki = 1:numel(clips)
        s = clips{ki}; N = numel(s);
        hs = cell(1, nF); lag = zeros(1, nF);
        for fi = 1:nF
            hs{fi} = sos_filter(F(fi).sos, s);
            lag(fi) = M.lag(s, hs{fi}, 160);
        end
        for seed = 1:NSEEDS
            t = t + 1; rng(1000*ci + 37*ki + seed);
            n = scale_noise(s, make_noise(cond{ci,2}, N, fs), cond{ci,3});
            x = s + n; snr0 = M.snr_db(s, x); mag0 = M.snr_mag(s, x);
            for fi = 1:nF
                yc = sos_filter(F(fi).sos, x);  yz = zero_phase(F(fi).sos, x);
                ys = {yc, M.align(yc, lag(fi)), yz};
                oref = {M.snr_oref(hs{fi}, sos_filter(F(fi).sos, n)), 0, ...
                        M.snr_oref(zero_phase(F(fi).sos, s), zero_phase(F(fi).sos, n))};
                oref{2} = oref{1};
                for mi = 1:nM
                    so = M.snr_db(s, ys{mi});
                    R.snr_out(ci,fi,mi,t) = so;  R.dsnr(ci,fi,mi,t) = so - snr0;
                    R.mse(ci,fi,mi,t) = M.mse(s, ys{mi});
                    R.oref(ci,fi,mi,t) = oref{mi} - snr0;
                    R.mag(ci,fi,mi,t) = M.snr_mag(s, ys{mi}) - mag0;
                end
            end
        end
    end
    fprintf('condition %d/%d done\n', ci, nC);
end
outdir = fullfile(here, '..', 'results'); if ~exist(outdir, 'dir'), mkdir(outdir); end
fid = fopen(fullfile(outdir, 'matlab_metrics_summary.csv'), 'w');
fprintf(fid, 'condition,filter,mode,snr_out,dsnr,dsnr_std,mse,dsnr_oref,dsnr_mag,n\n');
for ci = 1:nC, for fi = 1:nF, for mi = 1:nM
    d = squeeze(R.dsnr(ci,fi,mi,:));
    fprintf(fid, '%s,%s,%s,%.4f,%.4f,%.4f,%.6e,%.4f,%.4f,%d\n', cond{ci,1}, F(fi).name, modes{mi}, ...
        mean(squeeze(R.snr_out(ci,fi,mi,:))), mean(d), std(d), mean(squeeze(R.mse(ci,fi,mi,:))), ...
        mean(squeeze(R.oref(ci,fi,mi,:))), mean(squeeze(R.mag(ci,fi,mi,:))), nT);
end, end, end
fclose(fid); disp('wrote results/matlab_metrics_summary.csv');
