function F = design_filters(fs)
%DESIGN_FILTERS  The five filters compared in the paper (all as SOS matrices).
%   Returns a struct array with fields NAME and SOS.
if nargin < 1, fs = 16000; end
ny = fs/2;
F(1).name = 'LPF';     F(1).sos = butter_sos(4, 3400/ny, 'low');
F(2).name = 'HPF';     F(2).sos = butter_sos(4, 100/ny,  'high');
F(3).name = 'BPF';     F(3).sos = butter_sos(4, [100 3400]/ny, 'bandpass');
F(4).name = 'Notch';   F(4).sos = design_notch_sos(50, 30, fs);
F(5).name = 'Notch-H'; F(5).sos = design_notch_sos([50 100 150], 30, fs);
end
