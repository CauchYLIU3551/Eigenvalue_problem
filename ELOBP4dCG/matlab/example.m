load Kna2; load Mna2; load Hna2eig;
n = 1862; 

nb=4;
Nd=10;
m=3;
Z0=rand(2*n,nb);
shift=10;

E=eye(n,n);
E_plus=E;
E_minus=E;

opts.precond=1;
opts.precond_one=2;    % needed when opts.precond=1

[lamb,res,Z, hsty, info]=ELOBP4dCG(K, M, Z0, E, E, shift, Nd, m, 10^-5, opts);