% generate delay-time channel matrix (Eq. (4.55))
H_tilda=P*G*P.’
% generate delay-Doppler channel matrix (Eq. (6.1))
H=kron(Im,Fn)*(P’*G*P)*kron(Im,Fn’);