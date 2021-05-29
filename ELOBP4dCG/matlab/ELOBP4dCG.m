## Copyright (C) 2021 root
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} ELOBP4dCG (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: root <root@AFEPack01>
## Created: 2021-04-23

function [lambda, res, Z, hsty, info] = ELOBP4dCG (K, M, Z0, E_plus, E_minus, shift, Nd=10, m=3, tol=10^-3, opts)
%% This function computes the smallest k eigenvalues and corresponding eigenvec,
%% by ELOBP4dCG method.

%% Input:
%       K      array (n-by-n), SPD
%       M      array (n-by-n), SPD
%       E_plus      array (n-by-n), SPD
%       E_minus      array (n-by-n), SPD
%       Z0     array (2n-by-nb) whose columns span an approximate 
%              invariant subspace associated with (k0+1)st to (k0+k)th 
%              smallest positive eigenvalues
%       Nd,     int , the number of desired eigenvalues
%       m      int >=2, the order of krylov subspace default set as 3
%       tol    tolerence between true value and numerical result
%       opts   precondition options.
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

for i=1:nb
  RY(:,i)=CG(M,RY(:,i));
  RX(:,i)=CG(K,RX(:,i));
end

RY=[Y,RY];
RX=[X,RX];
Yw=M_GS(RY,M);
Xw=M_GS(RX,K);

[V1,S1,U1]=svd(Yw'*E_minus*Xw);

[s,idx]=sort(diag(S1),'descend'); 
ns=length(s); %s=s(1:nb);

X=Xw*U1(:,idx(1:nb)); nX=size(Xw,2); XmX=Xw(:,nb+1:nX)*U1(nb+1:nX,idx(1:nb));
Y=Yw*V1(:,idx(1:nb)); nY=size(Yw,2); YmY=Yw(:,nb+1:nY)*V1(nb+1:nY,idx(1:nb));
rho=ones(1,ns)./(s(1:ns))';
Rho=rho;

KX=K*X; 
MY=M*Y;

RX=KX-Y*diag(rho(1:nb)); RY=MY-X*diag(rho(1:nb));      % Residuals

nrmZ=sum(abs(X))+sum(abs(Y)); 
nrmR=(sum(abs(RX))+sum(abs(RY)))./(nrmH*nrmZ+rho(1:nb).*nrmZ);   % Residual errors in all-components

hsty.eig=rho(1:nb); hsty.res=nrmR
eig=rho(1:nb);
%R=nrmR;

% Convergence test
cvd_indx=find(nrmR<=tol); 
nk_cvd=length(cvd_indx); 
cvd_indx_c=find(nrmR>tol);

if nk_cvd > 0
   lam=[lam; (rho(cvd_indx))']; 
   res=[res nrmR(cvd_indx)]; 
   Z=[Z, [Y(:,cvd_indx); X(:,cvd_indx)]];
   V=[V, Y(:,cvd_indx)];
   n_cvd=n_cvd+nk_cvd;
   
   % make working block to nb again
   Xa=Xw*V1(:,idx(nb+1:nb+kk));           % Xa should have K-orthonormal columns
   Ya=Yw*U1(:,idx(nb+1:nb+kk));           % Ya should have M-orthonormal columns
   KXa=K*Xa+V*(shift*(V'*Xa));
   MYa=M*Ya;
   RXa=KXa-Ya*diag(rho(nb+1:nb+kk));
   RYa=MYa-Xa*diag(rho(nb+1:nb+kk));
   
   RX=[RX(:,cvd_indx_c) RXa]; RY=[RY(:,cvd_indx_c) RYa];
   rho=rho([cvd_indx_c nb+1:nb+kk]); 
   X=[X(:,cvd_indx_c) Xa]; Y=[Y(:,cvd_indx_c) Ya];
   KX=[KX(:,cvd_indx_c) KXa];               % Part of this KX does not include the updates for the ones that have just converged.
                                          % Probably doesn't matter because, they will be included in the
                                          % very next iteration and forward.
   MY=[MY(:,cvd_indx_c) MYa];
else
   rho=rho(1:nb);
end



iter=1;

while iter < max_iter && n_cvd<Nd
    XKX=X'*KX; % It should be I, but as precaution ...
    colXKX=chol(XKX); 
    invXKX=inv(colXKX); 
    X=X*invXKX;  % X'*shiftedK*X=I
 
    YMY=Y'*MY; % It should be I, but as precaution ...
    colYMY=chol(YMY); 
    invYMY=inv(colYMY); 
    Y=Y*invYMY;  % Y'*M*Y=I
    
    Xw=[]; Yw=[];
    for ie=1:m-1       

        % Preconditioning  
        if n_cvd == 0,            
           if opts.precond==1    
              if opts.precond_one==1
                 % P=K\RX, Q=M\RY
                 P=Rk\(Rk'\RX); Q=Rm\(Rm'\RY);
              elseif opts.precond_one==2
                 % P=K\RX, Q=M\RY
                 for i=1:nb,
                     %[P(:,i), error, iter, flag] = LinCG(K, CGx0, RX(:,i), CG_M, CGopts1); 
                     %[Q(:,i), error, iter, flag] = LinCG(M, CGx0, RY(:,i), CG_M, CGopts2); 
                     P(:,i)=CG(K, RX(:,i), CGx0);
                     Q(:,i)=CG(M, RY(:,i), CGx0);
                 end
              elseif opts.precond_one==3
                     % TBI
              end
           end               
        else % kc>0          
           if opts.precond==1
              if opts.precond_one==1
                 % P=[K+shift*V*V']\RX, Q=M\RY         
                 KinvRX=Rk\(Rk'\RX);
                 P=KinvRX-KiViT*(KiV'*RX);
                 Q=Rm\(Rm'\RY);
              elseif opts.precond_one==2
                 % P=[K+shift*V*V']\RX, Q=M\RY
                 for i=1:nb,
%                     [P(:,i), error, iter, flag] = LinCG(K, CGx0, RX(:,i), CG_M, CGopts1); 
%                     [Q(:,i), error, iter, flag] = LinCG(M, CGx0, RY(:,i), CG_M, CGopts2); 
                     P(:,i)=CG(K, RX(:,i), CGx0);
                     Q(:,i)=CG(M, RY(:,i), CGx0);
                 end
              elseif opts.precond_one==3
                    % TBI
              end
           end             
        end           
           
        Xw=[Xw P]; Yw=[Yw Q];
        if ie<m-1
           KX=K*P;    
           if n_cvd>0,
              KX=KX+V*(shift*(V'*P));
           end
           MY=M*Q;   
           RX=KX-Q*diag(rho); RY=MY-P*diag(rho);
        end
         
    end   
    
    Xn=MGSgS(K, shift, V, [Xw, XmX], X); Xw=[X, Xn];
    Yn=MGSg(M, [Yw, YmY], Y);            Yw=[Y, Yn];
    
    [U1,S1,V1] = svd(Yw'*Xw); [s,idx]=sort(diag(S1),'descend'); ns=length(s); %s=s(1:nb);
    X=Xw*V1(:,idx(1:nb)); nX=size(Xw,2); XmX=Xw(:,nb+1:nX)*V1(nb+1:nX,idx(1:nb));
    Y=Yw*U1(:,idx(1:nb)); nY=size(Yw,2); YmY=Yw(:,nb+1:nY)*U1(nb+1:nY,idx(1:nb));
    rho=ones(1,ns)./(s(1:ns))';  
 
    KX=K*X;  
    if n_cvd>0,
       KX=KX+V*(shift*(V'*X));
    end
    MY=M*Y;
    
    RX=KX-Y*diag(rho(1:nb)); RY=MY-X*diag(rho(1:nb));      % Residuals
    
    nrmZ=sum(abs(X))+sum(abs(Y)); 
    nrmR=(sum(abs(RX))+sum(abs(RY)))./(nrmH*nrmZ+rho(1:nb).*nrmZ);   % Residual errors in all-components
    hsty.eig=[hsty.eig; rho(1:nb)]; hsty.res=[hsty.res; nrmR];   
    
    % Convergence test
    cvd_indx=find(nrmR<=tol); nk_cvd=length(cvd_indx); 
    cvd_indx_c=find(nrmR>tol);
    % # of converged eigenpairs is kk.  Theorectically, all nb pairs could coverge at once. 
    % But that's unlikely, we assume they don't, i.e., kk<nb. 
    if nk_cvd>0
       lambda=[lambda; (rho(cvd_indx))']; res=[res nrmR(cvd_indx)]; 
       Z=[Z, [Y(:,cvd_indx); X(:,cvd_indx)]];
       V=[V, Y(:,cvd_indx)]; 
       n_cvd=n_cvd+nk_cvd;    
          
       if opts.precond==1
          if opts.precond_one==1
             % [K+shift*V*V']^{-1} = K^{-1}-shift*K^{-1}*V*[I+shift*V'*K^{-1}*V]^{-1}*V'*K^{-1}
             KiV=[KiV, Rk\(Rk'\Y(:,cvd_indx))];
             T=eye(n_cvd)+shift*V'*KiV; 
             invT=inv(T);
             KiViT=KiV*(shift*invT);
          elseif opts.precond_one==2
             CGopts1.update = 1;
             CGopts1.shift =shift;
             CGopts1.V = V;
          end
       end 
       
       % make working block to nb again
       Xa=Xw*V1(:,idx(nb+1:nb+nk_cvd));           % Xa should have K-orthonormal columns
       Ya=Yw*U1(:,idx(nb+1:nb+nk_cvd));           % Ya should have M-orthonormal columns
       KXa=K*Xa+V*(shift*(V'*Xa));
       MYa=M*Ya;
       RXa=KXa-Ya*diag(rho(nb+1:nb+nk_cvd));
       RYa=MYa-Xa*diag(rho(nb+1:nb+nk_cvd));
       
       RX=[RX(:,cvd_indx_c) RXa]; RY=[RY(:,cvd_indx_c) RYa];
       %rho=rho(kk+1:kk+nb); 
       rho=rho([cvd_indx_c nb+1:nb+nk_cvd]); 
       X=[X(:,cvd_indx_c) Xa]; Y=[Y(:,cvd_indx_c) Ya];
       KX=[KX(:,cvd_indx_c) KXa];               % Part of this KX does not include the updates for the ones that have just converged.
                                              % Probably doesn't matter because, they will be included in the
                                              % very next iteration and forward.
       MY=[MY(:,cvd_indx_c) MYa];
    else
       rho=rho(1:nb);
    end
    
    iter=iter+1; 
  
end   

info.itn = iter; 
Z=[Z [Y;X]];
endfunction
