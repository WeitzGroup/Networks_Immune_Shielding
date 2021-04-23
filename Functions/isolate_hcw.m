% This function takes the infected HCWs and moves them into an isolated
% class
% Function written by Adriana Lucia-Sanz
function new_node_status = isolate_hcw(node_status, is_hcw, pars)
% node_status (0 is S, 1 is E, 2 is I and 3 is R, and 4 infected but isolated)
 if ~(isequal(size(node_status), size(is_hcw)))
     error('Make sure the vectors node_status, is_pat, is_hcw are the same size')
 end
 
    new_node_status = node_status; %make a copy
    ind_I = find((node_status == 2) .* is_hcw); %index of infected HCW
    ind_I = ind_I(:).'; % make sure this is a row vector

    for i_I = ind_I
        if rand < pars.P_isolation
            new_node_status(i_I) = 4; % become isolated
        end
    end    
end