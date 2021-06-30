clc
clear
addpath Functions
addpath Data
addpath Figures
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

% Number of iterations
n_runs = 500;
R_init = [20,40,60,80,100]; % initial number of immunized indiv.
n_R_init = length(R_init);
% Initial conditions
inf_0 = 1; %initial number of infected hcws

% Vectors
Rinf_prew = zeros(n_runs, n_R_init);
Rinf_bas = zeros(n_runs, n_R_init);
%% BASELINE CASE 
for i_r = 1:n_R_init
    for i = 1:n_runs
        %Initial conditions with inf_0 randomly infected hcws and rec_0
        %randomly immunized individuals
        node_status = initial_cond(inf_0, R_init(i_r), pars);
        epi = zeros(4, N_steps);

        for i_t = 1:N_steps
            node_status = SEIR_stochastic_fct(adj, node_status, pars); %SEIR dynamics
            epi(:, i_t) = type_to_count(node_status); % Count each type
        end
        %save final recovered for run i and initial cond i_r
        Rinf_bas(i, i_r) = epi(4, end);
    end
end
%% PREWIRING CASE
for i_r = 1:n_R_init  
    for i = 1:n_runs
        node_status = initial_cond(inf_0, R_init(i_r), pars);      
        %prewire
        adj_new = pre_rewire(adj, node_status, is_pat, is_hcw);
        epi = zeros(4, N_steps);

        for i_t = 1:N_steps
            node_status = SEIR_stochastic_fct(adj_new, node_status, pars); %SEIR dynamics
            epi(:, i_t) = type_to_count(node_status); % Count each type
        end
        %save final recovered for run i and initial cond i_r
        Rinf_prew(i, i_r) = epi(4, end);
    end
end
%% Figure 4 prob. of outbreak 
for i_r = 1:n_R_init
    %matrix for prob outbreak for baseline and prewire
    P_outbreak = zeros(n - R_init(i_r), 2);
    for x = 0:(n - R_init(i_r) - inf_0)
        P_outbreak(x+1, 1) = sum((Rinf_bas(:,i_r) - R_init(i_r))>=x)./n_runs;
        P_outbreak(x+1, 2) = sum((Rinf_prew(:,i_r) - R_init(i_r))>=x)./n_runs;
    end
    
    figure(5 + i_r);
    hold on;
    plot(linspace(0,200, n - R_init(i_r)), P_outbreak(:, 1),...
        'Color',[0.5,0.5,0.5],'LineWidth',2.5)
    plot(linspace(0,200, n - R_init(i_r)), P_outbreak(:, 2),...
        'Color',[0.9,0.6,0.6],'LineWidth',2.5)  
    hold off;
    legend('no interventions','prewire','box','off')
    xlabel('Outbreak size');
    ylabel('Prob. outbreak');
    set(gca, 'FontSize',16);    
    box on
end
%% Figure 5 number of infections per run

for i_r = 1:n_R_init
    [a1,b1] = hist(Rinf_bas(:, i_r) - R_init(i_r) - inf_0, 0:3:200);
    [a2,b2] = hist(Rinf_prew(:, i_r) - R_init(i_r) - inf_0, 0:3:200);
    
    figure('Position', [200, 200, 750, 300]);
    hold all
   
    bar(b1, a1./n_runs,'LineWidth',1.5,...
        'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.7);
    bar(b2, a2./n_runs,'LineWidth',1.5,...
        'EdgeColor','none','FaceColor',[0.9,0.6,0.6],'FaceAlpha',0.7);
    hold off
    
    set(gca, 'YScale', 'log', 'FontSize',16);
    xlabel('Outbreak size');
    ylabel('Frequency');
    title(strcat(num2str(R_init(i_r)*100/n),'% Immunization'));
    box on
    legend('no interventions','prewire','box','off')
end