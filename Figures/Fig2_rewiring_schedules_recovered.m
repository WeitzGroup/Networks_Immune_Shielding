clc
clear
load epi_params

%Time parameters
T_vec = pars.T_vec;
N_steps = pars.N_steps;
n = pars.n; % network size
pars.n_A = 100;
load rand_bip_n200.mat; % load a bipartite network of 200 nodes called 'adj'
load hcw_pat_id.mat % load logical column vectors for hcw and pat
%Epi is the network with the number of S, E, I and R in every time step
zero_vector = pars.zero_vector;

%Rewire frequency
daily = 24*6;
weekly = daily*7;
rewire_freq = weekly;
n_runs = 500;
inf_0 = 10; %initial number of infected hcws
rec_0 = 0; %initial number of recovered individuals
%% No rewiring 
S_nah = zeros(n_runs, N_steps);
E_nah = zeros(n_runs, N_steps);
I_nah = zeros(n_runs, N_steps);
R_nah = zeros(n_runs, N_steps);

% Save cumulative of new infections per time per run
new_infections_baseline = zeros(N_steps, n_runs);

for i = 1:n_runs
    cumulative_inf=[];
    %Initial conditions
    node_status = initial_cond(inf_0, rec_0, pars);
    epi_temp = zeros(N_steps, 4);
    
    for i_t = 1:N_steps
        node_status = SEIR_stochastic_fct(adj, node_status, pars);
        epi_temp(i_t, :) = type_to_count(node_status);
        
        % calculate new infections
        inf_nodes = find(node_status == 2);
        cumulative_inf = [cumulative_inf; inf_nodes];
        cumulative_inf = unique(cumulative_inf);
        new_infections_baseline(i_t, i) = numel(cumulative_inf);
    end
    R_nah(i, :) = epi_temp(:, 4);
    I_nah(i,:) = epi_temp(:,3);
    S_nah(i, :) = epi_temp(:,1);
    E_nah(i,:) = epi_temp(:,2);
end
new_infections_baseline = new_infections_baseline';
new_infections_baseline = new_infections_baseline - inf_0;

%% How rewiring delays affect the dynamics of recovered nodes
freq_rewiring = [1 3 5 7]; % frequency of rewiring, daily, every 3 days, every 5 days, and weekly
R_all_rewfreq = zeros(n_runs, N_steps, length(freq_rewiring));
new_inf_freq = zeros(N_steps, n_runs, length(freq_rewiring));


for f = 1:length(freq_rewiring)
    rewire_freq = 24 * 6 * freq_rewiring(f);
    S_all = zeros(n_runs, N_steps);
    E_all = zeros(n_runs, N_steps);
    R_all = zeros(n_runs, N_steps);
    I_all = zeros(n_runs, N_steps);
   
    % Save cumulative of new infections per time per run
    new_infections_rew = zeros(N_steps, n_runs);
    
    for i = 1:n_runs 
        cumulative_inf=[];
        adj_new = adj; %for every run start with the same adj matrix
        node_status = initial_cond(inf_0, rec_0, pars);
        epi_temp = zeros(N_steps, 4);
        for i_t = 1:N_steps
            node_status = SEIR_stochastic_fct(adj_new, node_status, pars);
            epi_temp(i_t, :) = type_to_count(node_status);
            
            % calculate new infections
            inf_nodes = find(node_status == 2);
            cumulative_inf = [cumulative_inf; inf_nodes];
            cumulative_inf = unique(cumulative_inf);
            new_infections_rew(i_t, i) = numel(cumulative_inf);
            
            if mod(i_t, rewire_freq) == 0
                adj_new = rewire_all(adj_new, node_status, is_pat, is_hcw);
            end
        end
        R_all(i, :) = epi_temp(:, 4);
        I_all(i,:) = epi_temp(:,3);
        S_all(i, :) = epi_temp(:, 1);
        E_all(i,:) = epi_temp(:,2);
    end
    R_all_rewfreq(:, :, f) = R_all;
    new_inf_freq(:,:,f) = new_infections_rew;
end

new_inf_freq = new_inf_freq - inf_0;

%% Figure rewiring delays effects on outbreak size
steps = 1:(24*6):N_steps;
error_times = 10*(24*6):10*(24*6):N_steps;
gray_col = linspace(50,200,8);
col = repmat(gray_col', 1, 3);
col = col./255;
figure(4)
plot_rec(T_vec, n, new_infections_baseline, steps, error_times, col(8,:))
hold on
plot_rec(T_vec, n, new_inf_freq(:,:,4)', steps, error_times, col(7,:)) % weekly rewiring
plot_rec(T_vec, n, new_inf_freq(:,:,3)', steps, error_times, col(5,:)) % every 5 days rewiring
plot_rec(T_vec, n, new_inf_freq(:,:,2)', steps, error_times, col(3,:)) % every 3 days rewiring
plot_rec(T_vec, n, new_inf_freq(:,:,1)', steps, error_times, col(1,:)) % daily rewiring
hold off
legend('no rewiring', 'weekly', 'every 5 days', 'every 3 days', 'daily', 'location', 'best')
legend box off
title('Rewiring schedules and their effects on outbreak size')