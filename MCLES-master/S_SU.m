function [S, D] = S_SU(X, k, issymmetric)
if nargin < 3
    issymmetric = 1;
end;
if nargin < 2
    k = 5;
end;
[~, n] = size(X);
D = L2_distance_1(X, X);
[~, idx] = sort(D, 2); 
S = zeros(n);
for i = 1:n
    id = idx(i,2:k+2);
    di = D(i, id);
    S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;
if issymmetric == 1
    S = (S+S')/2;
end;