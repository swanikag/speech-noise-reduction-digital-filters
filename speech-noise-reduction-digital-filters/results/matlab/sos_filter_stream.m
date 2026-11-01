function [y, state] = sos_filter_stream(sos, xframe, state)
%SOS_FILTER_STREAM  Real-time (frame-by-frame) causal biquad cascade with carried state.
%   state is a K-by-2 matrix of filter states (initialise with zeros(K,2)).
y = xframe(:);
for k = 1:size(sos, 1)
    [y, state(k, :)] = filter_state(sos(k, 1:3), sos(k, 4:6), y, state(k, :));
end
end

function [y, zf] = filter_state(b, a, x, zi)
[y, zf] = filter(b, a, x, zi(:));
zf = zf(:).';
end
