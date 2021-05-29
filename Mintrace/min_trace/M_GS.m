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
## @deftypefn {} {@var{retval} =} M_GS (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: root <root@AFEPack01>
## Created: 2021-01-20

#################################################################################
## The function is used to compute the M-orthogonal modified Gram-Schmidt process
## Input: matrix B;
## Output: M-orthogonal matrix T;



function T = M_GS (B,M)
  [m,n]=size(B);
  T=zeros(m,n);
  T(:,1)=B(:,1);
  for i=2:n
    T(:,i)=B(:,i);
    for j=1:i-1
      T(:,i)=T(:,i)-projM(T(:,j),B(:,i),M);
    end
  end
  for k=1:n
    T(:,k)=T(:,k)/sqrt(T(:,k)'*M*T(:,k));
  end
  
endfunction
