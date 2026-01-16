% Method 1 (Eqs. (4.24) and (4.27))
Y_tilda=reshape(r,M,N);
Y=Y_tilda*Fn;
% Method 2 (Eq. (4.35))
y=kron(eye(M),Fn)*(P.’)*r;
Y=reshape(y,N,M).’;
% Method 3 (Eq. (4.35))
y=(P.’)*kron(Fn,eye(M))*r;
Y=reshape(y,N,M).’;