z=exp(1i*2*pi/N/M);
delay_spread=max(l_i);
% Generate discrete-time baseband channel in TDL form (Eq. (2.22))
gs=zeros(delay_spread+1,N*M);
for q=0:N*M-1
for i=1:taps
gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z∧(k_i(i)*(q-l_i(i)));
end
end
% Generate discrete-time baseband channel matrix (Eq. (4.38))
G=zeros(N*M,N*M);
for q=0:N*M-1
for ell=0:delay_spread
if(q>=ell)
G(q+1,q-ell+1)=gs(ell+1,q+1);
end
end
end