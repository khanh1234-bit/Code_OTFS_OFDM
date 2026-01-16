% Method 1: Using the TDL model (Eq. (4.36))
r=zeros(N*M,1);
for q=0:N*M-1
for ell=0:(delay_spread-1)
if(q>=ell)
r(q+1)=r(q+1)+gs(ell+1,q+1)*s(q-ell+1);
end
end
end
% Method 2: Using the time-domain channel matrix (G) (Eq. (4.37))
r=G*s.’;
% Method 3: Using the delay-time channel matrix (H_tilda) (Eq. (4.54))
x_tilda=reshape(X_tilda.’,N*M,1);
y_tilda=H_tilda*x_tilda;
r=P*y_tilda;
% Method 4: Using the delay-Doppler channel matrix (H) (Eq. (4.59))
x=reshape(X.’,N*M,1);
y=H*x;
r=P*kron(Im,Fn’)*y;