function [H] =  S_HT(M,P,S,F,V,alpha)
N = size(S,2);
I = eye(N);
AW = F*V.^(1/2);
A_syl = P'*P;
B_syl = alpha*AW*AW';
C_syl = -P'*M;
H = lyap(A_syl,B_syl,C_syl);
end