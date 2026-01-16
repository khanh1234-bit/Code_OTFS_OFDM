z=exp(1i*2*pi/N/M);
delay_spread=max(l_i);
l_cp=delay_spread;
% Generate discrete-time baseband channel in TDL form (Eq. (2.22))
gs=zeros(delay_spread+1,N*(M+l_cp));
for q=0:N*(M+l_cp)-1
for i=1:taps
gs(l_i(i)+1,q+1)=gs(l_i(i)+1,q+1)+g_i(i)*z∧(k_i(i)*(q-l_i(i)));
end
end
% Generate discrete-time baseband channel matrix (Eq. (4.93))
G_cp=zeros(N*M,N*M);
for n=0:N-1
for m=0:M-1
for ell=0:delay_spread
G_cp(m+n*M+1,n*M+mod(m-ell,M)+1)=gs(ell+1,m+n*M+l_cp+1);
end
end
end
% generate received signal after discarding CP per block
r=G_cp*s.’;