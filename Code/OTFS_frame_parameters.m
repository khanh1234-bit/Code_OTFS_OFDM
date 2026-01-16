% number of Doppler bins (time slots)
N=16;
% number of delay bins (subcarriers)
M=64;
% normalized DFT matrix
Fn=dftmtx(N);
Fn=Fn/norm(Fn);
% subcarrier spacing
delta_f=15e3;
% block duration
T=1/delta_f;
% carrier frequency
fc=4e9;
% speed of light
c=299792458;
% OTFS grid delay and Doppler resolution
delay_resolution = 1/(M*delta_f);
Doppler_resolution = 1/(N*T);