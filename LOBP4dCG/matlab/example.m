## get the matrix K and M from the data files.
data1=load('Kna2.mat');
data2=load('Mna2.mat');

# get the dimension of the problem;
K=data1.K;
M=data2.M;
n=max(size(K));

[eig_val, eig_vec,P,Q]=LOBP4dCG(K,M,4);