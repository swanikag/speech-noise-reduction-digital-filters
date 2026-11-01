function n = scale_noise(s, n, snr_db)
%SCALE_NOISE  Scale noise so that 10*log10(sum(s.^2)/sum(n.^2)) equals snr_db exactly.
n = n * sqrt(sum(s.^2) / (sum(n.^2) * 10^(snr_db/10)));
end
