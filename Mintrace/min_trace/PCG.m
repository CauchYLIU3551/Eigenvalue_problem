## Copyright (C) 2021 root
## 

## Author: cyliu 
## Created: 2021-01-20

# This function is used to solve the Rayleigh quetient minimization problem;
# The generalized problem is : Ax=lambda Mx,
# Output: minimal eigenvector x and the minimal eigenvalue;

function x = PCG (A, M, x0)
  tol=10^-5;
  n=size(M);
  x0=zeros(n,1);
  x0(1)=sqrt(1/M(1,1));
  x=x0;
  v=A*x;
  u=M*x;
  rho=v'*x/(u'*x);
  g0=2(v-rho*u);
  g1=g0;
  k=0;
  while norm(g)>tol
    if k==0
      p0=-g0;
      p1=p0;
      k=1
    else
      p0=p1;
      p1=-g1+g1'*M*g1/(g0'*M*g0)*p0;
    end
    
    ## Determine the smallest Ritz value rho and corresponding Ritz vector x of 
    ## (A,M) in R([x_k-1, p_k]);
    delta=find_step(A, M, x, p);
    x=x+delta*p;
    
    
    v=A*x;
    u=M*x;
    rho=x'*v/(x'*u);
    g0=g1;
    g1=2*(v-rho*u);
  end
    
end

