%{
 This is a generating function for different types of bipartite networks:
 nested (triangular)
 identity (I_n)
 complete (K_n)
 random
 regular with degree k
 Watts-Strogatz

Arguments are size of first partition (HCWs), second parition (Patients) 
and type of matrix 
Returns adjacency matrix adj and ID 
Written by Adriana Lucia-Sanz
Edited by Andreea Magalie
%}

function [adj,ID] = bipartite_nw(node1,node2,type)

% Nested network
if strcmp(type,'nested')
    bip=triu(ones(node1,node2),0); % triangular matrix 

% One-to-one network
elseif strcmp(type, 'one-to-one')
    bip=eye(node1,node2); % identity matrix

% All-to-all network
elseif strcmp(type, 'clique')
    bip=ones(node1,node2); %complete matrix

% Random network
elseif strcmp(type,'random')
    k=10; % node degree
    p=2*k/(node1+node2); % edge probability
    bip = (p>=rand(node1,node2)); % random bipartite graph with average degree k

% Regular network with degree k
elseif strcmp(type,'regular')
    k=10; % node degree
    bip=zeros(node1,node2);
    
    %This connects node i frome node1 with the next k nodes from node2
    for i=1:node1
        spand=i+k-1;
        if spand<=node2
            bip(i,i:spand)=1;
        %if node i is too large, wraps around 
        else
            spand=spand-node2;
            bip(i,i:end)=1;
            bip(i,1:spand)=1;
        end   
    end
elseif strcmp(type,'Watts-Strogatz')
    k = 10;
    P_connection_ws = 0.01;
    adj = WattsStrogatz(node1,k/2,P_connection_ws);
    print('Input type of network is invalid')
end

%created a (node1 + node2) x (node1 + node2) adj matrix out of bip 
if ~strcmp(type,'Watts-Strogatz')
    adj=[zeros(node1,node1),bip;bip',zeros(node2,node2)]; 
end

%creates node IDs
ID=cell(1,node1+node2);
for x=1:node1+node2
    ID(x)={strcat(['H_',num2str(x)])};
    if x>node1
        ID(x)={strcat(['P_',num2str(x-node1)])};
    end        
end

end
