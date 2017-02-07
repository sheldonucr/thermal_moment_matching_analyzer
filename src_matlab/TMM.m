function [ve] = TMM(timestamp, G, C, T0, power, ENV_);
warning off Matlab:nearlySingularMatrix;

[length,tnode]=size(power);                       % number of simulation points
time=timestamp(2:length,1)-timestamp(1:length-1,1); % time step = time(i+1)-time(i)

temp_seg = T0';
q=7;
ve(1,1:tnode)=T0;

[poles, pm] = mm(G, C, q,  ENV_, tnode);  
% B=zeros(tnode,1); B(3,1)=1;
% [poles,pm]=prima(G,diag(C),B,q);

for index=1:length-1
    [residue, moment0] = mm_2 (poles, pm, G, C, (power(index,:))', temp_seg, q, ENV_, tnode);
    
    seg_temp = zeros(tnode,1);
    for comp=1:tnode
        seg_temp(comp)=0;
        for ii=1:size(poles,1)
            seg_temp(comp)=seg_temp(comp)+residue(comp,ii)*exp(poles(ii)*time(index));
        end;
        seg_temp(comp) = moment0(comp) + seg_temp(comp) + ENV_; 
    end;
    ve(index+1,:)=seg_temp';
    temp_seg = seg_temp;
end;