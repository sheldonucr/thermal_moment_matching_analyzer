function [ residue, moment0] = mm_2(poles, pm, b, c_ver, power_mean, temp_seg, q, ENV_, TNode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute moments

%disp(poles);
q=size(poles,1);
G=b;
C = diag(c_ver);

%Am = inv(G)*C;
Am = G\C;
X0 = temp_seg-ENV_;
Bm = power_mean;

%moment0 = inv(G)*Bm;
moment0 = G\Bm;
moment(:,1) = moment0 - X0;
for i=2:q
    moment(:,i) = -Am*moment(:,i-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute residues
for component = 1:length(c_ver) 
    residue (component,:) = (pm \ moment(component,1:q)')';
end;

return;
