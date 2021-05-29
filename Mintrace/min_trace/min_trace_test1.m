% This is a test example for trace minimization problem;
% It will solve the generalized eigenvalue problem such as:
% Ax=lambda M x; 
% By the theorem of the trace minimization, we can get the p smallest eigenvalues
% of the matrix A corresponding to M. So this .m file will implement that.
% Input: matrix A and the related matrix M, the number of the eigenvalues to 
%        solve;
% Output: the summand of the p eigenvalues;

n=10;
p=3;
max_iter=20;
A=diag(2*ones(10,1));
A=A+diag(ones(9,1),-1)+diag(ones(9,1),1);
M=diag([1,2,3,4,5,6,7,8,9,10]);
V=zeros(n,p);
V(1,1)=1;
V(2,2)=sqrt(2)/2;
V(3,3)=sqrt(3)/3;
e=eye(p);

for k=1:max_iter
  W=A*V;  
  H=V'*W;
  [U,theta]=eig(H);
  X=V*U;
  R=W*U-M*X*theta;
  P=eye(n)-M*X/(X'*M*M*X)*X'*M;
  for i=1:p
    x=zeros(n,1);
    ## Maybe the PCG solver is used to solve the linear equations, not the 
    ## Rayleigh quotient minimization algorithm!
    #x=PCG(P*(A-theta(i,i)*M)*P,P*R*e(:,i),x);
    x=CG(P*(A-theta(i,i)*M)*P,P*R*e(:,i),x);
    X(:,i)=X(:,i)-x;
  end

  #V=M_GS(X,R);
  V=M_GS(X,M);
end