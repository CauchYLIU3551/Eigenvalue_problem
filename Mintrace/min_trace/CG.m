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

function x = CG (A, b, x0)
  x=x0;
  g=b-A*x;
  p=-g;
  test=A*p;
  delta=p'*g/(p'*A*p)
  x=x+delta*p;
  g=b-A*x;
  %num=1;
  while norm(g)>10^-5%&num<200
    beta=g'*A*p/(p'*A*p);
    p=-g+beta*p;
    delta=p'*g/(p'*A*p);
    x=x+delta*p;
    g=b-A*x;
    %num=num+1;
   end
   %num

endfunction
