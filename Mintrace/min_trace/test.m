n=10;
p=2;
max_iter=1;
tol=0.01;
A=diag(2*ones(10,1));
A=A+diag(ones(9,1),-1)+diag(ones(9,1),1);
M=diag([1,2,3,4,5,6,7,8,9,10]);
V=zeros(n,p);
for i=1:p
  V(i,i)=sqrt(i)/i;
end

e=eye(p);

for k=1:max_iter
  W=A*V;  
  H=V'*W;
  [U,theta]=eig(H);
  X=V*U;
  R=W*U-M*X*theta;
  if norm(R)<tol
    printf("the error satisfies the tolerence");
    break;
  end
  P=eye(n)-M*X/(X'*M*M*X)*X'*M;
  delta=zeros(n,p);
  for i=1:p
    x=zeros(n,1);
    ## Maybe the PCG solver is used to solve the linear equations, not the 
    ## Rayleigh quotient minimization algorithm!
    ####x=PCG(P*(A-theta(i,i)*M)*P,P*R*e(:,i),x);
    #x=CG(P*(A-theta(i,i)*M)*P,P*R*e(:,i),x);
    X*e(:,i)
    tempAX=A*X*e(:,i)
    tempPAX=P*A*X*e(:,i)
    %a=X'*M*ones(10,1)
   % b=(X'*M*M*X)
    %x
    %P*ones(10,1)
    

    x=CG(P*A*P,P*A*X*e(:,i),x);
    delta(:,i)=x;
    X(:,i)=X(:,i)-x;
  end
  #V=M_GS(X,R);
  V=M_GS(X,M);
  end

eigenvalue=diag(X'*A*X);