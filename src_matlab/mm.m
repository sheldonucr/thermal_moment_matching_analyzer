function [poles, pm] = mm(b, c_ver, q, ENV_, TNode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute poles
G = b;
C = diag(c_ver);

U1 = zeros(TNode, 1);
U1(3,1) = 1;

B = zeros(TNode, TNode);
B(3,3) = 1;

s_0 = B * U1;

U = zeros(TNode, (2*q));
U(:,1) = s_0;

[G2, C2, B2, M, V] = eks(q, C, G, B, U);

Eig = eig(G2^(-1)*C2);
 q = length(Eig);

 for i=1:q
     poles(i) = -1/Eig(i);
 end;
poles = poles';

for i = 1:q
    for j = 1:q
            pm(i,j) = -1 / (poles(j)^(i-1));
    end;
end;

return;