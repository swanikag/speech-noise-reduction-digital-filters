%PLOT_FIGURES  MATLAB/Octave versions of the main paper figures. Saves PNG (300 dpi) to ../figures_matlab
%   Run after run_experiments.m (needs ../results/matlab_metrics_summary.csv for the bar chart).
clear; close all;
if exist('OCTAVE_VERSION', 'builtin'), pkg load signal; end
here = fileparts(mfilename('fullpath')); addpath(here);
outdir = fullfile(here, '..', 'figures_matlab'); if ~exist(outdir, 'dir'), mkdir(outdir); end
fs = 16000; F = design_filters(fs); nF = numel(F);
col = [0 .45 .70; .84 .37 0; 0 .62 .45; .80 .47 .65; .90 .62 0];
sv = @(h, name) print(h, fullfile(outdir, [name '.png']), '-dpng', '-r300');

% ---------- Fig 1: filter responses (magnitude, causal vs zero-phase, group delay)
f = logspace(log10(20), log10(8000), 4000);
H = zeros(nF, numel(f));
for k = 1:nF
    H(k, :) = sos_resp(F(k).sos, f, fs);
end
h1 = figure('Position', [100 100 1100 330], 'visible', 'off');
subplot(1,3,1); hold on;
for k = 1:nF, semilogx(f, 20*log10(max(abs(H(k,:)), 1e-6)), 'Color', col(k,:), 'LineWidth', 1.2); end
set(gca, 'XScale', 'log'); grid on; axis([20 8000 -60 5]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend({F.name}, 'Location', 'south'); title('(a) Magnitude');
subplot(1,3,2); hold on;
for k = 1:2
    m = 20*log10(max(abs(H(k,:)), 1e-6));
    plot(f, m, '-', 'Color', col(k,:), 'LineWidth', 1.2); plot(f, 2*m, '--', 'Color', col(k,:), 'LineWidth', 1.2);
end
set(gca, 'XScale', 'log'); grid on; axis([20 8000 -60 5]); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend({'LPF causal','LPF zero-phase','HPF causal','HPF zero-phase'}, 'Location', 'south'); title('(b) Causal vs zero-phase');
fl = linspace(20, 4000, 40000); subplot(1,3,3); hold on;
for k = 1:nF
    h = sos_resp(F(k).sos, fl, fs);
    gd = -gradient(unwrap(angle(h)), 2*pi*fl) * 1000;
    plot(fl, gd, 'Color', col(k,:), 'LineWidth', 1.2);
end
set(gca, 'XScale', 'log'); grid on; axis([20 4000 0 8]); xlabel('Frequency (Hz)'); ylabel('Group delay (ms)'); title('(c) Group delay');
sv(h1, 'fig_filter_response');

% ---------- example trial: first clip, WGN 5 dB
[clips, names] = load_clips(fullfile(here, '..', 'data')); M = metrics();
s = clips{1}; N = numel(s); t = (0:N-1)'/fs;
rng(2026); n = scale_noise(s, make_noise('white', N, fs), 5); x = s + n;
y = cell(1, nF);
for k = 1:nF
    lag = M.lag(s, sos_filter(F(k).sos, s), 160);
    y{k} = M.align(sos_filter(F(k).sos, x), lag);
end

% ---------- Fig 2: waveforms
h2 = figure('Position', [100 100 560 650], 'visible', 'off');
sig = {s, x, y{1}, y{3}}; lab = {'Clean', 'Noisy (WGN, 5 dB)', 'LPF (causal)', 'BPF (causal)'};
for i = 1:4
    subplot(4,1,i); plot(t, sig{i}, 'k', 'LineWidth', 0.4); axis([0 3 -0.6 0.6]); grid on;
    ylabel('Amplitude'); text(2.95, 0.5, lab{i}, 'HorizontalAlignment', 'right');
    if i < 4, set(gca, 'XTickLabel', []); end
end
xlabel('Time (s)'); sv(h2, 'fig_waveforms');

% ---------- Fig 3: spectra (Welch PSD)
h3 = figure('Position', [100 100 1100 380], 'visible', 'off');
subplot(1,2,1); hold on;
sg = {s, x, y{1}, y{3}}; nm = {'Clean', 'Noisy', 'LPF', 'BPF'}; cc = {'k', [.55 .55 .55], col(1,:), col(3,:)};
for i = 1:4
    [p, fw] = welch_psd(sg{i}, 1024, fs); plot(fw, 10*log10(p + 1e-20), 'Color', cc{i}, 'LineWidth', 1);
end
grid on; xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)'); legend(nm, 'Location', 'southwest'); title('(a) 5 dB WGN'); xlim([0 8000]);
rng(7); nh = scale_noise(s, make_noise('hum', N, fs), 5); xh = s + nh;
subplot(1,2,2); hold on;
sh = {xh, sos_filter(F(4).sos, xh), sos_filter(F(5).sos, xh), sos_filter(F(2).sos, xh)};
nm2 = {'Noisy (hum, 5 dB)', 'Notch (50 Hz)', 'Notch-H', 'HPF (100 Hz)'}; c2 = {[.55 .55 .55], col(4,:), col(5,:), col(2,:)};
for i = 1:4
    [p, fw] = welch_psd(sh{i}, 4096, fs); m = fw <= 400; plot(fw(m), 10*log10(p(m) + 1e-20), 'Color', c2{i}, 'LineWidth', 1);
end
grid on; xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)'); legend(nm2, 'Location', 'southeast'); title('(b) 5 dB hum'); xlim([0 400]);
sv(h3, 'fig_spectra');

% ---------- Fig 4: spectrograms
h4 = figure('Position', [100 100 1100 330], 'visible', 'off');
sg = {s, x, y{3}}; tl = {'Clean', 'Noisy (WGN, 5 dB)', 'BPF output (causal)'}; vmax = [];
for i = 1:3
    subplot(1,3,i); [D, ff, tt] = spec_db(sg{i}, fs, 512, 128);
    if isempty(vmax), vmax = max(D(:)); end
    imagesc(tt, ff/1000, D, [vmax-70 vmax]); axis xy; colormap(jet); title(tl{i}); xlabel('Time (s)');
    if i == 1, ylabel('Frequency (kHz)'); end
end
sv(h4, 'fig_spectrograms');

% ---------- Fig 5: Delta-SNR bars (reads results CSV written by run_experiments.m)
csvf = fullfile(here, '..', 'results', 'matlab_metrics_summary.csv');
if exist(csvf, 'file')
    fid = fopen(csvf); C = textscan(fid, '%s %s %s %f %f %f %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 1); fclose(fid);
    conds = unique(C{1}, 'stable'); fn = {F.name}; D = zeros(numel(conds), nF); O = D;
    for ci = 1:numel(conds), for k = 1:nF
        idx = strcmp(C{1}, conds{ci}) & strcmp(C{2}, fn{k}) & strcmp(C{3}, 'causal_aligned');
        D(ci,k) = C{5}(idx); O(ci,k) = C{8}(idx);
    end, end
    h5 = figure('Position', [100 100 1100 400], 'visible', 'off');
    subplot(1,2,1); b = bar(D); for k = 1:nF, set(b(k), 'FaceColor', col(k,:)); end
    set(gca, 'XTick', 1:numel(conds), 'XTickLabel', conds); grid on; ylabel('\DeltaSNR, waveform (dB)'); title('(a) Causal, delay-compensated');
    subplot(1,2,2); b = bar(O); for k = 1:nF, set(b(k), 'FaceColor', col(k,:)); end
    set(gca, 'XTick', 1:numel(conds), 'XTickLabel', conds); grid on; ylabel('\DeltaSNR, output-referred (dB)'); legend(fn, 'Location', 'northwest');
    title('(b) Output-referred');
    sv(h5, 'fig_dsnr_causal');
else
    disp('Run run_experiments.m first to create the bar chart.');
end
disp(['figures written to ' outdir]);
