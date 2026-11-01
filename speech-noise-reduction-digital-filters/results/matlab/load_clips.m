function [clips, names] = load_clips(folder, target_rms)
%LOAD_CLIPS  Read all 16 kHz mono WAV files, remove DC, scale to RMS = target_rms.
if nargin < 2, target_rms = 0.1; end
d = dir(fullfile(folder, '*.wav')); names = {d.name}; clips = cell(1, numel(d));
for i = 1:numel(d)
    [x, fs] = audioread(fullfile(folder, d(i).name));
    assert(fs == 16000, 'expected 16 kHz');
    x = x(:, 1); x = x - mean(x);
    clips{i} = x * (target_rms / sqrt(mean(x.^2)));
end
end
