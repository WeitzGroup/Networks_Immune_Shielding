clc
clear
load epi_params

%Time parameters
T_vec = pars.T_vec;
N_steps = pars.N_steps;
n = pars.n; % network size
pars.n_A = 100;
% Network names/iterations
network=["rand_bip_n200.mat", "reg_bip_n200.mat", "ws_rand_n200.mat",...
    "ws_reg_n200.mat"];

load hcw_pat_id.mat % load logical column vectors for hcw and pat

%Epi is the network with the number of S, E, I and R in every time step
zero_vector = pars.zero_vector;

%Rewire frequency
daily = 24*6;
weekly = daily*7;

% Initial conditions
n_runs = 50;
inf_0 = 10;
rec_0 = 0; 

% Vectors
epi = zeros(5, N_steps); % 1 - S (0), 2 - E (1),3 - I(2), 4 - R(3), 5 - ISOLATED(4)

%% BASELINE CASE
Rinf_bas = zeros(n_runs, length(network));
for i = 1:length(network)
    for j = 1:n_runs
        load(network(i));
        adj_new = full(adj);  
        
        node_status = initial_cond(inf_0, rec_0, pars);        
        epi = zeros(4, N_steps);

        for i_t = 1:N_steps
            node_status = SEIR_stochastic_fct(adj_new, node_status, pars); %SEIR dynamics
            epi(:,i_t) = type_to_count(node_status); % Count each type
        end
    
        Rinf_bas(j, i) = epi(4,end);
    end
end

%% REWIRING CASE DAILY

Rinf_rew_daily = zeros(n_runs, length(network));

for i = 1:length(network)
    for j = 1:n_runs
        load(network(i));
        adj_new = full(adj);
        
        node_status = initial_cond(inf_0, rec_0, pars);
        epi = zeros(4, N_steps);
        
        for i_t = 1:N_steps    
            node_status = SEIR_stochastic_fct(adj_new, node_status, pars); %SEIR dynamics
            
            if rem(i_t, daily) == 0 % daily rewiring
                adj_new = rewire_all(adj_new, node_status, is_pat, is_hcw);
            end
            epi(:, i_t) = type_to_count(node_status);
        end
        Rinf_rew_daily(j, i)=epi(4, end);   
    end
end


%% REWIRING CASE WEEKLY
Rinf_rew_weekly = zeros(n_runs, length(network));

for i = 1:length(network)
    for j = 1:n_runs
        load(network(i));
        adj_new = full(adj);  
        
        node_status = initial_cond(inf_0, rec_0, pars);
        epi = zeros(4, N_steps);
        
        for i_t = 1:N_steps    
            node_status = SEIR_stochastic_fct(adj_new, node_status, pars); %SEIR dynamics
            
            if rem(i_t, weekly) == 0 % weekly rewiring
                adj_new = rewire_all(adj_new, node_status, is_pat, is_hcw);
            end
            epi(:, i_t) = type_to_count(node_status);
        end
        Rinf_rew_weekly(j, i)=epi(4, end);   
    end
end

%% Supplementary Figure 3a - Outbreak size for different network structures
% set(gcf, 'Position',  [200, 200, 1200, 450])% set position, width and height of plot
hold all;

data=[Rinf_bas(:,1),Rinf_rew_weekly(:,1),Rinf_rew_daily(:,1),Rinf_bas(:,2),...
    Rinf_rew_weekly(:,2),Rinf_rew_daily(:,2), Rinf_bas(:,3),Rinf_rew_weekly(:,3),...
    Rinf_rew_daily(:,3), Rinf_bas(:,4),Rinf_rew_weekly(:,4),Rinf_rew_daily(:,4)];

data=data-inf_0.*ones(size(data)); % remove initial infected

color_boxplot = brewermap(3, 'Set1');
boxplot(data,'color', color_boxplot)
ylabel('Outbreak size')
xlh = xlabel('Networks', 'Position', [4 -10]);
set(findobj(gca,'type','line'),'linew', 2)
set(gca,'XTickLabel',{' '})
text(1.2, -5, 'random', 'fontsize', 12, 'fontweight', 'bold')
text(3.1, -5, 'regular', 'fontsize', 12, 'fontweight', 'bold')
text(5.0, -5, 'ws-random', 'fontsize', 12, 'fontweight', 'bold')
text(8.9, -5, 'ws-regular', 'fontsize', 12, 'fontweight', 'bold')
set(gca, 'fontsize', 13, 'fontweight', 'bold')
ylim([0,200]);

