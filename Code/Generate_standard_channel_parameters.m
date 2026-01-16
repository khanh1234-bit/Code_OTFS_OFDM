% choose channel model, for example: EVA
delays=delays_EVA;
pdp=pdp_EVA;
% dB to linear scale
pdp_linear = 10.∧(pdp/10);
% normalization
pdp_linear = pdp_linear/sum(pdp_linear);
% number of propagation paths (taps)
taps=length(pdp);
% generate channel coefficients (Rayleigh fading)
g_i = sqrt(pdp_linear).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
% generate delay taps (assuming integer delay taps)
l_i=round(delays./delay_resolution);