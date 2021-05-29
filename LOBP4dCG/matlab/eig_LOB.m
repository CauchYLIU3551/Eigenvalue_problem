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
## @deftypefn {} {@var{retval} =} eig_LOB (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: root <root@AFEPack01>
## Created: 2021-04-19

function [eigenvector, eigenvalue] = eig_LOB (A,k)
  [orig_vec,orig_val]=eig(A);
  eig_val=diag(orig_val);
  num=max(size(eig_val));
  %temp_index=zeros(num/2);
  temp_index=[];
  
  % get the index array of the positive eigenvalues;
  for i=1:num
    if real(eig_val(i))>0
      temp_index(end+1)=i;
    endif
  end
  
  % get the k smallest positive eigenvalues index and using them to get
  % corresponding eigenvectors.
  eigenvalue=zeros(k,1);
  eigenvector=zeros(num,k);
  flag=max(eig_val);
  for i=1:k
    [val,ind]=min(eig_val(temp_index));
    eigenvalue(i)=val;
    eigenvector(:,i)=orig_vec(:,temp_index(ind));
    temp_index(ind)=[]; % once the entry is used then delete it;
  end
  
endfunction
