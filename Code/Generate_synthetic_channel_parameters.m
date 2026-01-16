% number of propagation paths
taps=6;
% maximum normalized delay and Doppler spread
l_max=4;
k_max=4;
% generate channel coefficients (Rayleigh fading) with uniform pdp
g_i = sqrt(1/taps).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
% generate delay taps uniformly from [0,l_max]
l_i = [randi([0,l_max],1,taps)];
l_i= l_i-min(l_i);
% generate Doppler taps (assuming uniform spectrum [-k_max,k_max])
k_i = k_max-2*k_max*rand(1,taps);