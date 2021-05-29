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

## Author: root <root@AFEPack01>
## Created: 2021-01-21

############################################################################
## This function is used to get the step size in the direction
## p_k. So that we can get the next iterate x_k+1=x_k+delta*p_k;

## Input: matrix A, related matrix M, iterate x, the direction p;

function delta = find_step (A, M, x, p)
  matrixA=[x'*A*x, x'*A*p; p'*A*x, p'*A*p];
  matrixM=[x'*M*x, x'*M*p; p'*M*x, p'*M*p];
  
  delta=eig(matrixA, matrixM)(1);
endfunction
