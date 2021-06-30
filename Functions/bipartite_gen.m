%This function generates a bipartite network with tot_E edges and n_A nodes
%in the first partition and n_B nodes in the second partition
function adj_new = bipartite_gen(n_A, n_B, tot_E)
N = n_A + n_B;
adj_new = zeros(N, N);
p = tot_E/(n_A * n_B);

if p >= 1
    error('n_A and n_B cannot form a bipartite network with 100 edges')
end

for i =1:n_A
    for j = 1:n_B
        if rand < p
            adj_new(i, n_A + j) = 1;
            adj_new(n_A + j, i) = 1;
        end
    end
end
end