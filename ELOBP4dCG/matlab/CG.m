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
## @deftypefn {} {@var{retval} =} CG (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: root <root@AFEPack01>
## Created: 2021-01-21

function x = CG (A, b, x0=zeros(max(size(A)),1),tol=10^-3)
  max_iter=max(size(A));
  x=x0;
  Ax=A*x;
  r=Ax-b;
  p=-r;
%  g=b-Ax;
%  p=-g
%  %test=A*p;
  Ap=A*p;
%  delta=p'*g/(p'*Ap)
%  x=x+delta*p
%  g=b-Ax
  num=0;
  while norm(r)>tol&&num<max_iter
    delta=-p'*r/(p'*Ap);
    x=x+delta*p;
    r=A*x-b;
    beta=r'*Ap/(p'*Ap);
    p=-r+beta*p;
    Ap=A*p;
%    delta=p'*g/(p'*Ap)
%    x=x+delta*p;
%    g=b-A*x;
    num=num+1;
   end
   %num

endfunction
