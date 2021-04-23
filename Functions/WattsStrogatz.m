function adj= WattsStrogatz(N,K,P_conn)
% Adj = WattsStrogatz(N,K,P_conn) returns a Watts-Strogatz model graph with
% N nodes, N*K edges, mean node degree 2*K, and rewiring probability P_conn.
%
% P_conn= 0 is a ring lattice, and P_conn = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.

s = repelem((1:N)',1,K); % creates a matrix of NxK entries which entries 
                         % for a row a_i=i
t = s + repmat(1:K,N,1); %repmat  gives a matrix NxK
                         % columns a_j=j to K. The sum gives a_ij=i+j
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta=P_conn
for source=1:N    
    switchEdge = rand(K, 1) < P_conn;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

G = graph(s,t);
adj=adjacency(G);
end