numerical_result4=[51.0402
49.9355
23.8171
11.0202
11.0871
9.71445e-16];
numerical_result3=[45.5041
45.1197
21.1185
10.1848
10.1972
2.22552e-16];
numerical_result2=[40.473
40.4536
19.9824
9.93224
9.9304
-7.11103e-15];
numerical_result1=[39.7268
39.7214
19.8008
9.88501
9.88497
3.93088e-15];
numerical_result5=[
39.5416
39.5401
19.7545
9.87341
9.87346
-2.88131e-14];
exact_result=[4*pi^2;4*pi^2;2*pi^2;pi^2;pi^2;0];
numerical_result=[numerical_result5,numerical_result1,numerical_result2,numerical_result3,numerical_result4];
h=[0.025,0.05,0.1,0.25,0.5];
error=numerical_result-exact_result;
order=zeros(6,4);
order2=zeros(1,4);
temp=max(error);
for i=1:4
  order(:,i)=log(error(:,i)./error(:,i+1))/log(h(i)/h(i+1));
  order2(i)=log(temp(i)/temp(i+1))/log(h(i)/h(i+1));
end
