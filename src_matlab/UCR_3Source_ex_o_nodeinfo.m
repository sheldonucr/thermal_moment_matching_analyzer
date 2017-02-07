


%Node(i).Type = 1;           %0 is interior or power node, 1 is prescribed ambient node, 2 is floating ambient node
%Node(i).ConnList = [2];     %Node connection list
%Node(i).R = [0.5351];       %Resistance values assoicated w/ node conn list
%Node(i).C = 0;              %Node capacitance
%Node(i).T0 = 35;            %Initial temperature of the node


NetInfo.NumNodes = 11;       %Number of nodes in the network

%Paste new network data below after optimization is completed
Node(1).Type = 1; 
Node(1).ConnList = [ 2]; 
Node(1).R = [ 10.039]; 
Node(1).C = 0; 
Node(1).T0 = 25; 

Node(2).Type = 0; 
Node(2).ConnList = [ 1 3]; 
Node(2).R = [ 10.039 11.73]; 
Node(2).C = 1.063987e-002; 
Node(2).T0 = 30; 

Node(3).Type = 0; 
Node(3).ConnList = [ 2 4]; 
Node(3).R = [ 11.73 0.010783]; 
Node(3).C = 6.812889e-003; 
Node(3).T0 = 40; 

Node(4).Type = 0; 
Node(4).ConnList = [ 3 5 8]; 
Node(4).R = [ 0.010783 1.5224 7.5708]; 
Node(4).C = 1.646101e+001; 
Node(4).T0 = 30; 

Node(5).Type = 0; 
Node(5).ConnList = [ 4 6]; 
Node(5).R = [ 1.5224 10.814]; 
Node(5).C = 1.726605e+001; 
Node(5).T0 = 40; 

Node(6).Type = 0; 
Node(6).ConnList = [ 5 7]; 
Node(6).R = [ 10.814 11.809]; 
Node(6).C = 2.728646e+000; 
Node(6).T0 = 30; 

Node(7).Type = 1; 
Node(7).ConnList = [ 6]; 
Node(7).R = [ 11.809]; 
Node(7).C = 0; 
Node(7).T0 = 25; 

Node(8).Type = 0; 
Node(8).ConnList = [ 4 9]; 
Node(8).R = [ 7.5708 0.00086476]; 
Node(8).C = 7.262033e+000; 
Node(8).T0 = 40; 

Node(9).Type = 0; 
Node(9).ConnList = [ 8 10]; 
Node(9).R = [ 0.00086476 0.010321]; 
Node(9).C = 1.147330e+002; 
Node(9).T0 = 35; 

Node(10).Type = 0; 
Node(10).ConnList = [ 9 11]; 
Node(10).R = [ 0.010321 3.2186]; 
Node(10).C = 5.061579e+001; 
Node(10).T0 = 25; 

Node(11).Type = 1; 
Node(11).ConnList = [ 10]; 
Node(11).R = [ 3.2186]; 
Node(11).C = 0; 
Node(11).T0 = 25; 

