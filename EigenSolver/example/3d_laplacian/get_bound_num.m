[m,n]=size(M);
N=0;
for i=1:m
  num = 0;
  for j=1:n
    if abs(M(i,j))>10^-6
      num = num+1;
    endif
  endfor
  if num==1
    N=N+1;
  end
end
