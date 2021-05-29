n=100;
K=zeros(n,n);
for i=1:n-1
  K(i+1,i)=-1;
  K(i,i+1)=-1;
  K(i,i)=2;
end
K(n,n)=2;
M=diag(1:n);

[eig_val, eig_vec,P,Q]=LOBP4dCG(K,M,4);
PP=[zeros(n,n),K;M,zeros(n,n)];