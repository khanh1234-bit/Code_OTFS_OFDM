% average QAM symbol energy
Es = mean(abs(qammod(0:mod_size-1,mod_size).∧ 2));
% SNR=Es/noise power
SNR_dB = 25;
SNR=10.∧(SNR_dB/10)
% noise power
sigma_w_2=Es/SNR
% generate Gaussian noise samples with variance=sigma_w_2
noise = sqrt(sigma_w_2/2)*(randn(N*M,1) + 1i*randn(N*M,1));
% add AWGN to the received signal
r=r+noise;