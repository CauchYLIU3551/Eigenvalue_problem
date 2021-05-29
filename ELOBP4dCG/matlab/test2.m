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


%[lamb,res,Z, rho]=ELOBP4dCG(K, M, Z0, E, E, shift, Nd, m, 10^-3 ,opts);


n=max(size(K));
nb=size(Z0,2);   % Most time nb is a little larger than Nd
lambda=[];
res=[];
V=[];
Z=[];
nrmKM=[norm(K,1); norm(M,1)];
nrmH=max(nrmKM);
n_cvd=0;

max_iter=round(min(0.1*n,200)*(Nd/nb));

if opts.precond==1
   if opts.precond_one==1
      Rk=chol(K); Rm=chol(M);
   elseif opts.precond_one==2
      CGopts1.nitn=10;
      CGopts1.tol=1e-2;
      CGopts1.met=2;
      CGopts2.nitn=10;
      CGopts2.tol=1e-2;
      CGopts2.met=2;
      CG_M=0;
      CGx0=zeros(n,1);
      P=zeros(n,nb); Q=zeros(n,nb);
   elseif opts.precond_one==3
      % TBI
   end
end

if opts.precond==1
  if opts.precond_one==1
     KiV=[];
  elseif opts.precond_one==2
     CGopts1.update = 0;
     CGopts2.update = 0;
  elseif opts.precond_one==3
         % TBI
  end
end 

%initial approximate
X=Z0(n+1:2*n,:); 
Y=Z0(1:n,:);  

KX=K*X;
XKX=X'*KX; colXKX=chol(XKX); invXKX=inv(colXKX); X=X*invXKX; KX=KX*invXKX; % X'*shiftedK*X=I
MY=M*Y;
YMY=Y'*MY; colYMY=chol(YMY); invYMY=inv(colYMY); Y=Y*invYMY; MY=MY*invYMY; % Y'*M*Y=I

W=X'*Y; W=0.5*(W+W');
RX=KX*W-E_plus*Y; RY=MY*W-E_minus*X;

%for i=1:nb
%  RY(:,i)=CG(M,RY(:,i));
%  RX(:,i)=CG(K,RX(:,i));
%end
%
%RY=[Y,RY];
%RX=[X,RX];
%Yw=M_GS(RY,M);
%Xw=M_GS(RX,K);
 for i=1:nb,
     [P(:,i), error, iter, flag] = LinCG(K, CGx0, RX(:,i), CG_M, CGopts1); 
     [Q(:,i), error, iter, flag] = LinCG(M, CGx0, RY(:,i), CG_M, CGopts2); 
 end
 % expand X; note ||X||_K=1
Xw=MGSgS(K, shift, V, P, X); Xw=[X, Xw];

% expand Y; note ||Y||_M=1
Yw=MGSg(M, Q, Y); Yw=[Y, Yw];
[V1,S1,U1]=svd(Yw'*E_minus*Xw);

[s,idx]=sort(diag(S1),'descend'); 
ns=length(s); %s=s(1:nb);

X=Xw*U1(:,idx(1:nb)); nX=size(Xw,2); XmX=Xw(:,nb+1:nX)*U1(nb+1:nX,idx(1:nb));
Y=Yw*V1(:,idx(1:nb)); nY=size(Yw,2); YmY=Yw(:,nb+1:nY)*V1(nb+1:nY,idx(1:nb));
rho=ones(1,ns)./(s(1:ns))';