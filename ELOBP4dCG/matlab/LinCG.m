function [x, error, iter, flag] = LinCG(A, x, b, M, cgopts)
%
% LinCG.m solves the symmetric positive definite linear system Ax=b 
% using the Conjugate Gradient method with preconditioning.
%
% input   A        REAL symmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix, not used if cgopts.met  =  0, 2
%         cgopts   struct, options
%                  cgopts.nitn =  INTEGER maximum number of iterations
%                  cgopts.tol  =  REAL error tolerance
%                  cgopts.met  =  0, no preconditioner
%                              =  1, preconditioner = inv(M'*M)
%                              =  2, diagonal preconditioner
%
% output  x        REAL, solution vector
%         error    REAL, history for normalized residuals 
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it
%---------------------------------------------------------------

  flag = 0;                                 % initialization
  iter = 0;
  max_it=cgopts.nitn;
  tol=cgopts.tol;

  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), 
      bnrm2 = 1.0; 
      error = 0.0;
      x=zeros(size(x));
      return
  end

  r = b - A*x;
  error = norm( r ) / bnrm2;
  if ( error <= tol ) 
     return;
  end
  if cgopts.met == 2,
     M=diag(A);
  end

  for iter = 1:max_it                       % begin iteration

     if cgopts.met == 0
        z = r;
     elseif cgopts.met == 1,
        z  = M\( (M')\ r);
     elseif cgopts.met == 2,
        z  = r./M;
     end
     
     rho = (r'*z);

     if ( iter > 1 ),                       % direction vector
        beta = rho / rho_1;
        p = z + beta*p;
     else
        p = z;
     end

     q = A*p;
     alpha = rho / (p'*q );
     x = x + alpha * p;                    % update approximation vector

     r = r - alpha*q;                      % compute residual
     error = norm( r ) / bnrm2;            % check convergence
     if ( error <= tol ), break, end 

     rho_1 = rho;

  end

  if ( error > tol ) flag = 1; end         % no convergence

% END cg.m
