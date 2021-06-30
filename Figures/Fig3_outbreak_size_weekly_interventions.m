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
isolation_freq = weekly;
isolation_rewiring_freq = weekly;
n_runs = 500;
inf_0 = 10; %initial number of infected hcws
rec_0 = 0; %initial number of recovered individuals
%% No rewiring (Baeline case)
S_nah = zeros(n_runs, N_steps);
E_nah = zeros(n_runs, N_steps);
I_nah = zeros(n_runs, N_steps);
R_nah = zeros(n_runs, N_steps);

% Given that every run is stochastic (random) then we need to run the
% SEIR model multiple times to obtain an average behavior of the model. In this case, we
% run our model 200 times and obtain how the number of S, E, I, and R nodes
% changes through time

for i = 1:n_runs   
    cumulative_inf=[];
    node_status = initial_cond(inf_0, rec_0, pars);
    epi_temp = zeros(N_steps, 4);
    
    for i_t = 1:N_steps
        node_status = SEIR_stochastic_fct(adj, node_status, pars); % Run the SEIR model
        epi_temp(i_t, :) = type_to_count(node_status);
    end
    R_nah(i, :) = epi_temp(:, 4);
    I_nah(i,:) = epi_temp(:,3);
    S_nah(i, :) = epi_temp(:,1);
    E_nah(i,:) = epi_temp(:,2);
end

%% Rewiring intervention
S_all = zeros(n_runs, N_steps);
E_all = zeros(n_runs, N_steps);
R_all = zeros(n_runs, N_steps);
I_all = zeros(n_runs,N_steps);

% save node status by the last time point
for i = 1:n_runs
    cumulative_inf=[];
    adj_new = adj; %for every run start with the same adj matrix
    node_status = initial_cond(inf_0, rec_0, pars);    
    epi_temp = zeros(N_steps, 4);
    
    for i_t = 1:N_steps
        node_status = SEIR_stochastic_fct(adj_new, node_status, pars);
        epi_temp(i_t, :) = type_to_count(node_status);
   
        % Rewiring at 'x' frequency
        if mod(i_t, rewire_freq) == 0
            adj_new = rewire_all(adj_new, node_status, is_pat, is_hcw); % rewiring intervention
        end
    end
    R_all(i, :) = epi_temp(:, 4);
    I_all(i,:) = epi_temp(:,3);
    S_all(i, :) = epi_temp(:, 1);
    E_all(i,:) = epi_temp(:,2);
end

%% Isolation of HCWs
S_iso = zeros(n_runs, N_steps);
E_iso = zeros(n_runs, N_steps);
R_iso = zeros(n_runs, N_steps);
I_iso = zeros(n_runs,N_steps);


for i = 1:n_runs 
    cumulative_inf=[];
    node_status = initial_cond(inf_0, rec_0, pars);
    
    for i_t = 1:N_steps 
        node_status = SEIR_stochastic_fct(adj, node_status, pars);

        if mod(i_t, isolation_freq) == 0
             node_status = isolate_hcw(node_status, is_hcw, pars); % isolation intervention
        end
        S_iso(i, i_t) = sum(node_status == 0);
        E_iso(i, i_t) = sum(node_status == 1);
        I_iso(i, i_t) = sum(node_status == 2);
        R_iso(i, i_t) = sum(node_status == 3);  
    end
end

%% Isolation + rewiring

S_iso_rew = zeros(n_runs, N_steps);
E_iso_rew = zeros(n_runs, N_steps);
R_iso_rew = zeros(n_runs, N_steps);
I_iso_rew = zeros(n_runs,N_steps);


for i = 1:n_runs 
    cumulative_inf=[];
    adj_new = adj; %for every run start with the same adj matrix
    node_status = initial_cond(inf_0, rec_0, pars);
    
    for i_t = 1:N_steps
        node_status = SEIR_stochastic_fct(adj_new, node_status, pars);

        if mod(i_t, isolation_rewiring_freq) == 0
             node_status = isolate_hcw(node_status, is_hcw, pars); % first do isolation of infected hcws
             adj_new = rewire_all(adj_new, node_status, is_pat, is_hcw); % then do rewiring
        end
        S_iso_rew(i, i_t) = sum(node_status == 0);
        E_iso_rew(i, i_t) = sum(node_status == 1);
        I_iso_rew(i, i_t) = sum(node_status == 2);
        R_iso_rew(i, i_t) = sum(node_status == 3);  
    end
end

%% Figure: Boxplot comparing the outbreak size with different weekly interventions

% The outbreak size or epidemic size is the total number of people that got
% infected during an epidemic. That translates to the total of recovered people
% by the end of the simulation, because we have saved the number of
% recovered nodes through time for each of the four main interventions,
% then we can easily calculate the epidemic or outbreak size associated to
% each intervention.

% Number of recovered individuals at the end of each intervention
% Baseline (500 simulations)
%rec_baseline = new_infections_baseline(:,end)-inf_0;
rec_baseline = R_nah(:, end)-inf_0;


% weekly Isolation (500 simulations)
%rec_isolation = new_infections_iso(:,end)-inf_0;
rec_isolation = R_iso(:,end)-inf_0;


% weekly rewiring (500 simulations)
% rec_rewiring = new_infections_rew(:,end)-inf_0;
rec_rewiring = R_all(:, end)-inf_0;


% weekly Isolation + rewiring (500 simulations)
%rec_isorew = new_infections_isorew(:,end)-inf_0;
rec_isorew = R_iso_rew(:,end)-inf_0;


% create a 500 by 4 matrix containing the epidemic sizes for the above 4
% interventions, each row represent a single run (and we have 500 total
% runs)
data = [rec_baseline rec_isolation rec_rewiring rec_isorew];

figure(5)
color_boxplot = brewermap(5, 'Set1');
color_boxplot(1,:) = [0.7843    0.7843    0.7843]; % gray for baseline
color_boxplot(2,:) = [0.3020    0.6863    0.2902]; % green for isolation
color_boxplot(3,:) = [0.2157    0.4941    0.7216]; % blue for rewiring
boxplot(data,'color', color_boxplot)
ylabel('Outbreak size')
xlh = xlabel('Interventions', 'Position', [2.5   -12]);
set(findobj(gca,'type','line'),'linew', 2)
set(gca,'XTickLabel',{' '})
text(0.8, -5, 'Baseline', 'fontsize', 12, 'fontweight', 'bold')
text(1.8, -5, 'Isolation', 'fontsize', 12, 'fontweight', 'bold')
text(2.8, -5, 'Rewiring', 'fontsize', 12, 'fontweight', 'bold')
text(3.5, -5, 'Isolation + Rewiring', 'fontsize', 12, 'fontweight', 'bold')
set(gca, 'fontsize', 13, 'fontweight', 'bold')
title('Outbreak size for weekly interventions')

%% Figure: Probability of having an outbreak of size >= 'x'
% Calculate the probability of having an outbreak of size 'x' as the number
% of simulations with a final outbreak size greater than or equal to 'x'
% divided by the total number of simulations

outbreak_sizes = 1:200;
figure(6)
color_boxplot = brewermap(5, 'Set1');
color_boxplot(1,:) = [0.7843    0.7843    0.7843]; % gray for baseline
color_boxplot(2,:) = [0.3020    0.6863    0.2902]; % green for isolation
color_boxplot(3,:) = [0.2157    0.4941    0.7216]; % blue for rewiring
probability_outbreak = [sum(rec_baseline >= outbreak_sizes)/n_runs; sum(rec_isolation >= outbreak_sizes)/n_runs;...
    sum(rec_rewiring >= outbreak_sizes)/n_runs; sum(rec_isorew >= outbreak_sizes)/n_runs ];

for i = 1:size(probability_outbreak, 1)
    plot(outbreak_sizes, probability_outbreak(i,:), 'LineWidth', 2.5, 'color', color_boxplot(i,:))
    hold on
end
hold off
xlim([outbreak_sizes(1) outbreak_sizes(end)])
title('Probability of having an outbreak of size \geq x')
ylabel('Probability of outbreak')
xlabel('Outbreak size')
legend ('Baseline', 'Isolation', 'Rewiring', 'Isolation + Rewiring', 'location', 'best')
legend boxoff
set(gca, 'fontsize', 14)