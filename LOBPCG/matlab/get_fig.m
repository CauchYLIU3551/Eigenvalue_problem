n_cov = 5;

[eig, vec, residual] = LOBPCG (A, M, n_cov, tol=10^-6, max_iter=100);
 [x,y] = size(residual);

 
 for i = 1 : n_cov 
 plot(1:y,residual(i,:));
 hold on 
 endfor
 set(gca,'Yscale','log');
 xlabel('Iteration number');
 ylabel('Residual');
 hold on;
 
 labels = {};
 for i = 1 : n_cov
   switch i
   case 1
   labels(end+1) = sprintf("%dst eigenvalue", i);
   
   case 2
   labels(end+1) = sprintf("%dnd eigenvalue", i);
   
   case 3
   labels(end+1) = sprintf("%drd eigenvalue", i);
   
   otherwise
   labels(end+1) = sprintf("%dth eigenvalue", i);
   end
 endfor
 legend(labels);
 title("LOBPCG computing 5 eigenpairs");