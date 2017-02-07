function [G2, C2, B2, M, V] =  eks(q, C, G, B, U)
n = size(C, 1);
M = zeros(n, q);

M(:,1) = G \ (B * U(:, 1));
for i = 2: q
    M(:, i) = G \ (B * U(:, i) - C * M(:, i - 1));
end

Alpha = zeros(1, q);
R = zeros(n, q);
etol = 1e-15;

M(:, 1) = G \ (B * U(:, 1));
Alpha(1) = norm(M(:, 1));
R(:, 1) = M(:, 1) / Alpha(1);
alp = Alpha(1);

R(:, 2) = G \ (alp * B * U(:, 2) - C * R(:, 1));
h  = R(:, 1)' * R(:, 2);
R(:, 2) = R(:, 2) - h * R(:, 1);
Alpha(2) = norm(R(:, 2));
R(:, 2) = R(:, 2) / Alpha(2);
alp = alp * Alpha(2);

for i = 3 : q
    r_pre = zeros(n, 1);
    for j = 1 : i - 2
        h = R(:, j)' * R(:, i - 1);
        r_pre = r_pre + h * R(:, j);
    end;
    R(:, i) = G \ (alp * B * U(:, i) - C *(R(:, i-1) + Alpha(i-1) * r_pre));
    
    r_pre = zeros(n, 1);
    for j = 1: i - 1
        h  = R(:, j)' * R(:, i);
        r_pre = r_pre + h * R(:, j);
    end;
    R(:, i) = R(:, i) - r_pre;
    
    if(norm(R(:, i)) < etol)
        break;
    else
        Alpha(i) = norm(R(:, i));
        R(:, i) = R(:, i) / Alpha(i);
        alp = alp * Alpha(i);
    end;
end;
%--------------------------------------------------------------------------

if norm(R(:, i)) < 1
    q = i - 1;
end;

V = R(:, 1:q);
%--------------------------------------------------------------------------

G2 = V' * G * V;
C2 = V' * C * V;
B2 = V' * B;
