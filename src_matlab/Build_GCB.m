function [G,C,T0,np,nv,P_pos,V_pos] =  Build_GCB(NetInfo, Node)

NumOfNodes=NetInfo.NumNodes;
R=zeros(NumOfNodes,NumOfNodes);
nv=0; 
np=0;
for i=1:NumOfNodes          % parsing the netlist
    for j=1:size(Node(i).ConnList,2);
        R(i,Node(i).ConnList(j))=Node(i).R(j);
    end;
    if(Node(i).Type==0)     %0 is interior or power node, 1 is prescribed ambient node, 2 is floating ambient node
        np=np+1;
        P_pos(np,1)=i;
        T0(np)=Node(i).T0;
        C(np)=Node(i).C;
    else
        if(Node(i).Type==1)
            nv=nv+1; 
            V_pos(nv,1)=i; 
            V_pos(nv,2)=Node(i).ConnList(1);
            V_pos(nv,3)=Node(i).R(1);
        end;
    end;
end;

G=zeros(np,np);         %generate G matrix
k=1;
for x=1:np
    for y=1:np
        i=P_pos(x); j=P_pos(y);
        if(i==j)
            for k=1:NumOfNodes
                if(R(i,k)~=0)
                    G(x,y)=G(x,y)+1/R(i,k);
                end;
            end;
        else
            if(R(i,j)~=0) 
                G(x,y)=-1/R(i,j); 
            end;
        end;               
    end;
end;
            
% for i=1:nv
%     x=V_pos(i,2);
%     for j=1:np
%         if(P_pos(j)==x)
%             B(i)=j;
%             break;
%         end;
%     end;
% end;

