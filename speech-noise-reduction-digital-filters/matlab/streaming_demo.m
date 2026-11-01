%STREAMING_DEMO  Real-time (frame-by-frame) causal filtering; verifies equality with offline-batch causal output.
clear; if exist('OCTAVE_VERSION', 'builtin'), pkg load signal; end
here = fileparts(mfilename('fullpath')); addpath(here);
fs = 16000; [clips, names] = load_clips(fullfile(here, '..', 'data')); s = clips{1};
rng(2026); x = s + scale_noise(s, make_noise('white', numel(s), fs), 5);
F = design_filters(fs);
for fi = 1:numel(F)
    sos = F(fi).sos; batch = sos_filter(sos, x);
    for frame = [80 160 320 512]
        state = zeros(size(sos, 1), 2); y = zeros(0, 1); tt = [];
        for i = 1:frame:(numel(x) - frame + 1)
            tic; [yf, state] = sos_filter_stream(sos, x(i:i+frame-1), state); tt(end+1) = toc; %#ok<AGROW>
            y = [y; yf]; %#ok<AGROW>
        end
        fprintf('%-8s frame %3d (%4.1f ms): max|stream-batch| = %.2e, mean %.3f ms, max %.3f ms\n', ...
            F(fi).name, frame, 1000*frame/fs, max(abs(y - batch(1:numel(y)))), 1000*mean(tt), 1000*max(tt));
    end
end
