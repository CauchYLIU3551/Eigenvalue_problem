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
## @deftypefn {} {@var{retval} =} LOBP4dCG (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: root <root@AFEPack01>
## Created: 2021-04-19

function [eigenvalue, eigenvector,P,Q] = LOBP4dCG (K,M,k=1,max_iter=1000, eps=10^-5, PHI=0)
  n=max(size(K));
  % Compute the initial X_0 Y_0,
  X0=rand(n,k);
  Y0=rand(n,k);  
  Z=[Y0;X0];
  X1=X0;
  Y1=Y0;
  
  % build the precondition matrix PHI
  if PHI==0
    PHI=[zeros(n,n),eye(n,n);eye(n,n),zeros(n,n)];
  endif
  
  % begin iteration
  rho=zeros(k,1);
  i=0;
  normerror=1;
  while i<max_iter&& normerror>eps
    if i==0
      for j=1:k
        rho(j)=Rho(X0(:,j),Y0(:,j),K,M);
      endfor    
    endif
    P=K*X1-Y1*diag(rho);
    Q=M*Y1-X1*diag(rho);
    
    %[Q;P]=PHI*[P;Q];
    if i==0
      U=[X0,P];
      V=[Y0,Q];
    else
      U=[X1,X0,P];
      V=[Y1,Y0,Q];
    endif
    U=gshmit(U);
    V=gshmit(V);
    W=U'*V;
    W1=W';
    W2=eye(size(W));
    
    %% Build matrix H_SR
    H_SR=[zeros(size(W)),W1'^-1*U'*K*U*W1^-1;W2'^-1*V'*M*V*W2^-1,zeros(size(W))];
    size(H_SR);
    %% Compute the eigenpairs of H_SR;
    [temp_eig_vec,eigenvalue]=eig_LOB(H_SR,k);
   
    l=max(size(W2));
    
    X0=X1;Y0=Y1;
    x=temp_eig_vec((l+1):(2*l),:);
    y=temp_eig_vec(1:l,:);
    %% Update X1 and Y1
    X1=U*W1^-1*x;
    Y1=V*W2^-1*y;
    rho=eigenvalue;
    
    %% normalize X1 and Y1;
    %for j=1:k
      %X1(:,j)=X1(:,j)/norm(X1(:,j));
   %   Y1(:,j)=Y1(:,j)/norm(Y1(:,j));
    %endfor
    i=i+1;
    normerror=max(norm(P),norm(Q));
  endwhile  
  disp(["After ", num2str(i), " iterations, the algorithm converge."]);
  eigenvector=[Y1;X1];
endfunction
