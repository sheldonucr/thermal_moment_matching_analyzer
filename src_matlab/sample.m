clear;
clear all;

UCR_3Source_ex_o_nodeinfo;                     % netlist infor, with two variables Node, and NetInfo

[G,C,T0,np,nv,P_pos,V_pos]=Build_GCB(NetInfo, Node);        % get G,C and T0 matrix, and positin matrix.
    

load data.txt;                
[timestamp, sources]=Get_trace(data, np, nv, P_pos, V_pos);

ENV_=0;                                                     % global ENV_ (GND) is 0 degree.
[ve] = TMM(timestamp, G, C, T0, sources, ENV_);

plot(timestamp, ve);