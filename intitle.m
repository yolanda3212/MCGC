function [F,temp,lambda,idx, distX] =  intitle(X)

k = 15;
r = -1;
num = size(X,2);


distX = L2_distance_1(X,X);
%distX = sqrt(distX);
[distX1, idx] = sort(distX,2);
A = zeros(num);
rr = zeros(num,1);
for i = 1:num
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;

if r <= 0
    r = mean(rr);
end;
lambda = mean(rr);

A0 = (A+A')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;
[F, temp, evs]=eig1(L0);