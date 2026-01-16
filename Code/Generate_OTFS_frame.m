%modulation size
mod_size=4;
% number of information symbols in one frame
N_syms_per_frame=N*M;
% number of information bits in one frame
N_bits_per_frame=N*M*log2(mod_size);
% generate random bits
tx_info_bits=randi([0,1],N_bits_per_frame,1);
% QAM modulation
tx_info_symbols=qammod(tx_info_bits,mod_size,’gray’,’InputType’,’bit’);
% Generate the MxN OTFS delay-Doppler frame
X=reshape(tx_info_symbols,M,N);
% Vectorized OTFS frame information symbols
x=reshape(X.’,N*M,1);