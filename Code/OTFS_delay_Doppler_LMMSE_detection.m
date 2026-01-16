% vectorize Y
y=reshape(Y.’,N*M,1);
% estimated delay-Doppler matrix (Eq. (6.18))
x_hat=(H’*H+sigma_w_2)∧(-1)*(H’*y);
% QAM demodulation
x_hat=qamdemod(x_hat,mod_size,’gray’);