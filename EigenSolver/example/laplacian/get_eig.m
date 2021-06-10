% Compute the eigenpairs of the matrix A and M which are produced by AFEPack.
%This function is used to confirm the numerical results

filename = "stiff_matrix.txt";
filename2 = "mass_matrix.txt";

delimiterIn = ' ';

A = importdata(filename,delimiterIn);
M = importdata(filename2, delimiterIn);

[u,v] = eig(A,M);