z=exp(1i*2*pi/N/M);
delay_spread=max(l_i);
l_zp=delay_spread;
% Generate discrete-time baseband channel in TDL form (Eq. (2.22))
gs=zeros(delay_spread+1,N*(M+l_zp));
for q=0:N*(M+l_zp)-1
for i=1:taps
gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z∧(k_i(i)*(q-l_i(i)));
end
end
% Generate discrete-time baseband channel matrix (Eq. (4.109))
G_zp=zeros(N*M,N*M);
for n=0:N-1
    for m=0:M-1
for ell=0:delay_spread
if(m>=ell)
G_zp(m+n*M+1,m+n*M-ell+1)=gs(ell+1,m+n*M+l_zp+1);
end
end
end
end
% generate received signal after discarding ZP per block
r=G_zp*s.’;