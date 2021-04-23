% This functions reconnects R-R links in a network consisting of only 
% R and S individuals
% Function written by Adriana Lucia-Sanz adapted from rewire_all

function adj_new = pre_rewire(adj_now, node_status, is_pat, is_hcw)
%code will not function if these 3 vectors are not the same size
if ~(isequal(size(node_status), size(is_pat)) && isequal(size(node_status), size(is_hcw)))
    error('Make sure the vectors node_status, is_pat, is_hcw are the same size')
end
 
is_sus = (node_status == 0);
is_rec = (node_status == 3);
adj_new = adj_now;
%________________________________________________________
%__________Minimize R(Pat) - R(HCW) links________________
%__________by R(Pat) - S(HCW) links______________________
%________________________________________________________
rec_pat = find(is_rec.*is_pat); % recovered patients
sus_pat = find(is_sus.*is_pat); % susceptible patients

%if we have rec pat AND sus pat
if ~(isempty(rec_pat) || isempty(sus_pat))
    N_Rp = size(rec_pat, 1);
    N_Sp = size(sus_pat, 1);
    N_pair = N_Rp * N_Sp;
    for i_pair = 1:N_pair % Do pair extractions for N_pair times
        i_Rp = randi(N_Rp);
        i_Sp = randi(N_Sp);
        node_Rp = rec_pat(i_Rp);
        node_Sp = sus_pat(i_Sp);
        
        % Find all R HCWs linked to the R patient but not the S patient
        rec_hcw = find(is_rec.* is_hcw.* adj_new(:, node_Rp).* (1-adj_new(:, node_Sp)));
        % Find all S linked to the S patient but not the R patient
        sus_hcw = find((is_sus).* adj_new(:, node_Sp).* (1 - adj_new(:, node_Rp)));
        
        N_rec_h = size(rec_hcw, 1);
        N_sus_h = size(sus_hcw, 1);
        N_swap = min(N_rec_h, N_sus_h);
        
        % Perform the swaps
        if N_swap > 0
            % Generate a list of unique random numbers to select the swapped elements
            nodes_rec_h = rec_hcw(randperm(N_rec_h,N_swap));
            nodes_sus_h = sus_hcw(randperm(N_sus_h,N_swap));
            
            % Remove existing edges
            adj_new(node_Rp, nodes_rec_h)=0;
            adj_new(nodes_rec_h, node_Rp)=0;
            adj_new(node_Sp, nodes_sus_h)=0;
            adj_new(nodes_sus_h, node_Sp)=0;
            % Add new swapped edges
            adj_new(node_Sp, nodes_rec_h)=1;
            adj_new(nodes_rec_h, node_Sp)=1;
            adj_new(node_Rp, nodes_sus_h)=1;
            adj_new(nodes_sus_h, node_Rp)=1;
        end
        
    end
end