%Takes initial number of infected hcws inf_0, initial number of recovered/
%immunized individuals and pars

function node_status = initial_cond(inf_0, rec_0, pars)
node_status = pars.zero_vector;
init_R = randsample(1:pars.n, rec_0);
node_status(init_R) = 3;% recovered individuals
not_rec = find(node_status == 0);
hcws_id_not_rec = intersect(1:pars.n/2, not_rec);
init_I = randsample(hcws_id_not_rec, inf_0);
node_status(init_I) = 2; %infected HCW
end