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
## @deftypefn {} {@var{retval} =} gshmit (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: root <root@AFEPack01>
## Created: 2021-04-19

function [Q,R]=gshmit(A)
  % Input: V is an m by n matrix of full rank m<=n
  % Output: an m-by-n upper triangular matrix R
  % and an m-by-m unitary matrix Q so that A = Q*R.
  [m,n]=size(A);R=zeros(n);
  Q=zeros(m,n);
  R=zeros(n,n);
  for j=1:n
    v=A(:,j);
    for i=1:j-1
      R(i,j)=Q(:,i)'*A(:,j);
      v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
  endfor
end