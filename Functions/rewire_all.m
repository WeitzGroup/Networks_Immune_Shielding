%{ 
Rewiring algorithm - changes adj matrix of a network to minimizes number of
SI links while maintaining node degrees constant 
Input parameters: adj matrix, node status, is_pat, is_hcw
Returns: New adjacency matrix after rewiring
Function written by Joey Leung and Andreea Magalie inspired by Strona et al
(Nature comm 2014)
%}

function adj_new = rewire_all(adj_now, node_status, is_pat, is_hcw)
%code will not function if these 3 vectors are not the same size
 if ~(isequal(size(node_status), size(is_pat)) && isequal(size(node_status), size(is_hcw)))
     error('Make sure the vectors node_status, is_pat, is_hcw are the same size')
 end

%exposed considered susceptible 
is_sus = (node_status == 0) + (node_status == 1);
is_inf = (node_status == 2);
is_rec = (node_status == 3);

adj_new = adj_now;
%________________________________________________________
%________First we minimize I(Pat) - S(HCW) links_________
%________________________________________________________
inf_pat = find(is_inf.*is_pat); % Infected patients
rs_pat = [find(is_rec.*is_pat); find(is_sus.*is_pat)]; %R or S patients
%if we have inf pat AND rs pat
if ~(isempty(inf_pat) || isempty(rs_pat))
    N_Ip = size(inf_pat, 1);
    N_RSp = size(rs_pat, 1);
    N_pair = N_Ip * N_RSp;
    for i_pair = 1:N_pair % Do pair extractions for N_pair times
        i_Ip = randi(N_Ip); %select random I pat
        i_RSp = randi(N_RSp); %select random RS pat
        node_Ip = inf_pat(i_Ip);
        node_RSp = rs_pat(i_RSp);
        % Find all S HCWs linked to the I patient but not the RS patient
        sus_hcw = find(is_sus .* is_hcw .* adj_new(:, node_Ip) .* (1-adj_new(:, node_RSp)));
        % Find all R or I HCWs linked to the RS patient but not the I patient
        ri_hcw = find((is_rec+is_inf) .* adj_new(:, node_RSp) .* (1 - adj_new(:, node_Ip)));
        
        N_sus_h = size(sus_hcw, 1);
        N_ri_h = size(ri_hcw, 1);
        N_swap = min(N_sus_h, N_ri_h);
        
        % Perform the swaps
        if N_swap > 0
            % Generate a list of unique random numbers to select the swapped elements
            nodes_sus_h = sus_hcw(randperm(N_sus_h,N_swap));
            nodes_rec_h = ri_hcw(randperm(N_ri_h,N_swap));
            
            % Remove existing edges
            adj_new(node_Ip, nodes_sus_h)=0;
            adj_new(nodes_sus_h, node_Ip)=0;
            adj_new(node_RSp, nodes_rec_h)=0;
            adj_new(nodes_rec_h, node_RSp)=0;
            % Add new swapped edges
            adj_new(node_RSp, nodes_sus_h)=1;
            adj_new(nodes_sus_h, node_RSp)=1;
            adj_new(node_Ip, nodes_rec_h)=1;
            adj_new(nodes_rec_h, node_Ip)=1;
        end
        
    end
end

%________________________________________________________
%________Then we minimize I(HCW) - S(PAT) links__________
%________________________________________________________
inf_hcw = find(is_inf.*is_hcw); % Infected patients
rs_hcw = [find(is_rec.*is_hcw); find(is_sus.*is_hcw)]; %R and S patients
%If there are both infected and (R or S) hcws
if ~(isempty(inf_hcw)|| isempty(rs_hcw))
    N_Ihcw = size(inf_hcw, 1);
    N_RShcw = size(rs_hcw, 1);
    N_pair = N_Ihcw * N_RShcw;
    for i_pair = 1:N_pair % Do pair extractions for N_pair times
        i_Ihcw = randi(N_Ihcw);
        i_RShcw = randi(N_RShcw);
        node_Ihcw = inf_hcw(i_Ihcw);
        node_RShcw = rs_hcw(i_RShcw);
        
        % Find all S Patients linked to the I hcw but not the RS hcw
        sus_pat = find(is_sus .* is_pat .* adj_new(:, node_Ihcw) .* (1-adj_new(:, node_RShcw)));
        % Find all R or I patients linked to the RS hcw but not the I hcw
        ri_pat = find((is_rec+is_inf) .* adj_new(:, node_RShcw) .* (1 - adj_new(:, node_Ihcw)));
        
        N_sus_p = size(sus_pat, 1);
        N_ri_p = size(ri_pat, 1);
        N_swap = min(N_sus_p, N_ri_p);
        
        % Perform the swaps
        if N_swap > 0
            % Generate a list of unique random numbers to select the swapped elements
            nodes_sus_p = sus_pat(randperm(N_sus_p,N_swap));
            nodes_rec_p = ri_pat(randperm(N_ri_p,N_swap));
            
            % Remove existing edges
            adj_new(node_Ihcw, nodes_sus_p)=0;
            adj_new(nodes_sus_p, node_Ihcw)=0;
            adj_new(node_RShcw, nodes_rec_p)=0;
            adj_new(nodes_rec_p, node_RShcw)=0;
            % Add new swapped edges
            adj_new(node_RShcw, nodes_sus_p)=1;
            adj_new(nodes_sus_p, node_RShcw)=1;
            adj_new(node_Ihcw, nodes_rec_p)=1;
            adj_new(nodes_rec_p, node_Ihcw)=1;
        end
        
    end
end
end