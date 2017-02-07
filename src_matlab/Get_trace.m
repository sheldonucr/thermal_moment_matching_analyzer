function [timestamp, power]=Get_trace(data,np,nv,P_pos,V_pos)

length=size(data,1);
timestamp=data(:,1);
power=zeros(length,np);

for i=1:nv
    x=V_pos(i,1); 
    y=V_pos(i,2);
    data(:,y+1)=data(:,y+1)+data(:,x+1)/V_pos(i,3);
end;

for i=1:np
    power(:,i)=data(:,P_pos(i)+1);
end;
