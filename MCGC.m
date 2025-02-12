function [result] = MCGC(X, alpha,lambda, d, ~, maxIters, gt)
C = size(unique(gt),1);  
V = size(X,2);    
N = size(X{1},2); 
MAXiter = 1000;   
REPlic = 20;     
pn = 15;          
islocal = 1;
gamm =5; 
for i=1:V
    X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1); 
end
for i=1:V    
    D{i} = size(X{i},1);
end
SD = 0;
M = [];
for i=1:V
    SD = SD + D{i};    
    M = [M;X{i}];  
end 
G = cell(1,V);
HT = cell(1,V);
S = cell(1,V);
for i = 1:V
    [S{i}, ~] = S_SU(X{i}, pn, 0);
end
U = zeros(N);
for i = 1:V
    U = U + S{i};
end
U = U/V;
for j = 1:N
    U(j,:) = U(j,:)/sum(U(j,:));
end
w = ones(1,V)/V;
[F,temp] = intitle(M);
for it=1:maxIters 
    for v = 1 : V
        D{v} = size(X{v},1);
        G{v} = zeros(D{v},d);                           
        HT{v} = rand(d,N);       
        G{v} = UpdateG(HT{v},X{v},G{v});      
    end
    for v = 1:V
        HT{v} = S_HT(X{v},G{v},S{v},F,temp,alpha);
    end
    for v = 1:V 
       [S{v},~,gamma]= S_S1(HT{v},U,w(v),C,gamm/alpha);
    end
    for v = 1:V
        US = U - S{v};
        distUS = norm(US, 'fro')^2;
        if distUS == 0
            distUS = eps;
        end
        w(v) = 0.5/sqrt(distUS);
    end  
    dist = L2_distance_1(F',F');
    U = zeros(N);
    for i=1:N
        idx = zeros();
        for v = 1:V
            s0 = S{v}(i,:);
            idx = [idx,find(s0>0)];
        end
        idxs = unique(idx(2:end));
        if islocal == 1
            idxs0 = idxs;
        else
            idxs0 = 1:N;
        end
        for v = 1:V
            s1 = S{v}(i,:);
            si = s1(idxs0);
            di = dist(i,idxs0);
            mw = V*w(v);
            lmw = lambda/mw;
            q(v,:) = si-0.5*lmw*di;
        end
        U(i,idxs0) = S_U(q,V);
        clear q;
    end
    xC= U;
    xC= (xC+xC')/2;
    D1 = diag(sum(xC));
    L = D1-xC;
    [F, temp, ev]=eig1(L, C, 0);
  
   obj1=0;   obj2=0;   
   obj3=0;   obj4=0;   obj5=0;
   for v = 1:V
        obj1 = obj1+norm((X{v}-G{v}*HT{v}),'fro')^2;  
        obj2 = obj2+norm((HT{v}*F*temp.^(1/2)),'fro')^2; 
        obj3 = norm(S{v},'fro')^2;
        obj4 = lambda*trace(F'*L*F);
        obj5 = obj5 + w(v)*norm((U-S{v}),'fro')^2;
        Obj(it) = obj1+alpha*obj2+gamm*obj3+obj4+obj5; 
        if (it>1 && (abs(Obj(it)-Obj(it-1))/Obj(it-1)) < 10^-5)
            break;
        end
    
    end
    
end
l = kmeans(F,C,'maxiter',MAXiter,'replicates',REPlic);
%[fm, Precision, Recall] = compute_f(gt, l); 
result = ClusteringMeasure(gt,l);
