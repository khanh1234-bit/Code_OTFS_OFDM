% estimated time domain samples (Eq. (6.19))
s_hat=(G’*G+sigma_w_2)∧(-1)*(G’*r);
% MxN estimated delay-Doppler symbols (using Method 1 in code 12)
X_hat_tilda=reshape(s_hat,M,N);
X_hat=X_hat_tilda*Fn;
x_hat=reshape(X_hat.’,N*M,1);
% QAM demodulation
x_hat=qamdemod(x_hat,mod_size,’gray’);