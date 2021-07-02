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
## @deftypefn {} {@var{retval} =} LOBPCG (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: root <root@AFEPack01>
## Created: 2021-07-01

function [eigenvalues, eigenvectors] = LOBPCG (A, B, p, tol=10^-5, max_iter=100)
  [m,n] = size(A);
  x = rand(m,p);
  mu = zeros(p,1);
  r = zeros(m,p);
  x0 = zeros(m,p);
  for i=1:1%max_iter
    for j=1:p
      tempx = x(:,j);
      Ax = A*tempx;
      Bx = B*tempx;
      mu(j) = tempx'*Bx/(tempx'*Ax);
      r(:,j) = Bx - mu(j)*Ax;
    endfor
    r
    mu
    if max(max(abs(r)))<tol
      break
    endif
    if i == 1
      tempM = [r,x];
    elseif
      tempM = [r, x, x0];
    endif
    tempM=M_GS(tempM,A);
    [u,v] = eig(tempM'*B*tempM);
    [v,index] = sort(diag(v),'descend');
    x0 = x;
    x = tempM*u(:,index(1:p));  
  endfor
  eigenvalues = mu;
  eigenvectors = x;

endfunction
