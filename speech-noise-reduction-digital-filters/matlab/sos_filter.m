function y = sos_filter(sos, x)
%SOS_FILTER  Causal filtering with a cascade of biquads (zero initial state).
y = x(:);
for k = 1:size(sos, 1)
    y = filter(sos(k, 1:3), sos(k, 4:6), y);
end
end
