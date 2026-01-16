Im=eye(M);
% row-column permutation matrix (Eq. (4.33))
P=zeros(N*M,N*M);
for j=1:N
for i=1:M
E=zeros(M,N);
E(i,j)=1;
P((j-1)*M+1:j*M,(i-1)*N+1:i*N)=E;
end
end
% Method 1 (Eqs. (4.19) and (4.20))
X_tilda=X*Fn’;
s=reshape(X_tilda,1,N*M);
% Method 2 (Eq. (4.35))
s=P*kron(Im,Fn’)*x;
% Method 3 (Eq. (4.35))
s=kron(Fn’,Im)*P*x;