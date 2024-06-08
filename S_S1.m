% min_{A>=0, A*1=1, F'*F=I}  trace(D'*A) + r*||A||^2 + 2*lambda*trace(F'*L*F)
% written by Feiping Nie on 2/9/2014
function [A,beta,lambda] = S_S1(X,U,w,c,alpha, ~, r, islocal)
k = 15;
NITER = 30;
num = size(X,2);
if nargin < 7
    islocal = 1;
end;
if nargin < 6
    r = -1;
end;

distX = L2_distance_1(X,X);
[~, idx] = sort(distX,2);
A = zeros(num);
rr = zeros(num,1);
for i = 1:num
    id = idx(i,2:k+2);
    di = distX(i, id);
    numerator = di(k+1)-di+2*w*U(i,id(:))-2*w*U(i,id(k+1));   
    rr(i) = k*di(k+1)-sum(di(1:k));
    q = rr(i);
    denominator2 = 2*w*sum(U(i,id(1:k)))-2*k*w*U(i,id(k+1));
    A(i,id) = max(numerator/(rr(i)+denominator2+eps),0);    
    beta = (q-2*k*w*U(i,id(k+1))-2*w)/2;
end;

if r <= 0
    r = mean(rr);
end;
lambda = mean(rr);

A0 = (A+A')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;
[F, temp, evs]=eig1(L0, c, 0);

if sum(evs(1:c+1)) < 0.00000000001
    error('The original graph has more than %d connected component', c);
end;

for iter = 1:NITER
    distf = L2_distance_1(F',F');
    A3 = zeros(num);
    for i=1:num
        if islocal == 1
            idxa0 = idx(i,2:k+1);
        else
            idxa0 = 1:num;
        end;
        dfi = distf(i,idxa0);
        dxi = distX(i,idxa0);
        ad = -(2*alpha*dxi+lambda*dfi)/(2*r);
        A3(i,idxa0) = EProjSimplex_new(ad);
    end;

    A1 = (A3+A3')/2;
    D = diag(sum(A1));
    L = D-A1;
    F_old = F;
    [F, temp, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;

    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > 0.00000000001
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;  F = F_old;
    else
        break;
    end;

end;

[~, y]=conncomp(A1); 
y = y';

