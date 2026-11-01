function h = sos_resp(sos, f, fs)
%SOS_RESP  Complex frequency response (row vector) of a biquad cascade at frequencies f (Hz).
f = f(:).'; zi = exp(-1j*2*pi*f/fs); h = ones(size(f));
for r = 1:size(sos, 1)
    h = h .* (sos(r,1) + sos(r,2)*zi + sos(r,3)*zi.^2) ./ (sos(r,4) + sos(r,5)*zi + sos(r,6)*zi.^2);
end
end
