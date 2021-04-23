% This function takes adjacency matrix adj, node types (S (0), E (1), I (2), R(3) and isolated (4) in
% and transition probabilities p to change the state of every node in the network in one
% time step. 
function new_node_types = SEIR_stochastic_fct(adj, node_types, pars)
    new_node_types = node_types;

    ind_E = find(node_types == 1); % exposed nodes
    ind_I = find(node_types == 2); % infected nodes
    ind_HI = find(node_types == 2 | node_types == 4); %isolated and infected
    
    %make sure these are row vectors
    ind_E = ind_E(:).';
    ind_I = ind_I(:).';
    ind_HI = ind_HI(:).';

    %first change E nodes to I
    for i = ind_E
        if rand < pars.P_EI
            new_node_types(i) = 2;
        end
    end
    
    %second change isolated and infected nodes to R
    for i = ind_HI
        if rand < pars.P_IR
            new_node_types(i) = 3;
        end
    end
    
    %now look for potential I - S links
    for i = ind_I
        % if node i_I makes a contact
        if rand < pars.P_contact
            %chose a neighbor at random
            i_neigh = find(adj(i, :) == 1);
            i_inf = i_neigh(randi([1 length(i_neigh)]));
            %if link is to  a S, change S to E
            if node_types(i_inf) == 0
                new_node_types(i_inf) = 1;
            end
        end
    end
end