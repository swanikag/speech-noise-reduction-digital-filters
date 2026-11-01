function y = zero_phase(sos, x)
%ZERO_PHASE  Offline (non-causal) forward-backward filtering, no edge padding.
%   Needs the whole signal in memory and the FUTURE samples -> cannot run in real time.
y = sos_filter(sos, x(:));
y = flipud(sos_filter(sos, flipud(y)));
end
