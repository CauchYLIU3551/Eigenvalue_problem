% Compute the eigenpairs of the matrix A and M which are produced by AFEPack.
%This function is used to confirm the numerical results

filename = "stiff_matrix.txt";
filename2 = "mass_matrix.txt";
filename3 = "stiff_matrix2.txt";
filename4 = "mass_matrix2.txt";
delimiterIn = ' ';

A = importdata(filename,delimiterIn);
M = importdata(filename2, delimiterIn);
A2 = importdata(filename3,delimiterIn);
M2 = importdata(filename4,delimiterIn);

[u,v] = eig(A,M);