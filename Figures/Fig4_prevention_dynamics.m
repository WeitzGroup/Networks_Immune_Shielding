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

load rand_bip_n200.mat; % load a bipartite network of 200 nodes called 'adj'
load hcw_pat_id.mat % load logical column vectors for hcw and pat

%Epi is the network with the number of S, E, I and R in every time step
zero_vector = pars.zero_vector;
pars.is_hcw=is_hcw;

%Rewire frequency
daily = 24*6;
weekly = daily*7;
rewire_freq = weekly;

% Number of iterations
n_runs = 500;
pars.n_A = 100;
% Initial conditions
inf_0 = 1; %initial number of infected hcws
rec_0 = 60; %initial number of recovered individuals: 30% of the LTC

epi_S = zeros(n_runs, N_steps);
epi_E = zeros(n_runs, N_steps);
epi_I = zeros(n_runs, N_steps);
epi_R = zeros(n_runs, N_steps);


%% BASELINE CASE 
 
for i_iter = 1:n_runs
    %Initial conditions with inf_0 randomly infected hcws and rec_0
    %randomly immunized individuals
    node_status = initial_cond(inf_0, rec_0, pars);
    epi = zeros(4, N_steps);

    for i_t = 1:N_steps
        node_status = SEIR_stochastic_fct(adj, node_status, pars); %SEIR dynamics
        epi(:,i_t) = type_to_count(node_status); % Count each type
    end

    % Save each run 
    epi_S(i_iter,:) = epi(1,:);
    epi_E(i_iter,:) = epi(2,:);
    epi_I(i_iter,:) = epi(3,:);
    epi_R(i_iter,:) = epi(4,:); 

end

S_bas_mean = mean(epi_S);
S_bas_std  = std(epi_S);
E_bas_mean = mean(epi_E);
E_bas_std = std(epi_E);
I_bas_mean = mean(epi_I);
I_bas_std = std(epi_I);
R_bas_mean  = mean(epi_R);
R_bas_std = std(epi_R);


%% PREWIRING CASE
  
for i_iter = 1:n_runs
    node_status = initial_cond(inf_0, rec_0, pars);      
    %prewire
    adj_new = pre_rewire(adj, node_status, is_pat, is_hcw);
    epi = zeros(4, N_steps);

    for i_t = 1:N_steps
        node_status = SEIR_stochastic_fct(adj_new, node_status, pars); %SEIR dynamics
        epi(:,i_t) = type_to_count(node_status); % Count each type
    end

    % Save each run 
    epi_S(i_iter ,:) = epi(1,:);
    epi_E(i_iter ,:) = epi(2,:);
    epi_I(i_iter ,:) = epi(3,:);
    epi_R(i_iter ,:) = epi(4,:);

end

S_prew_mean = mean(epi_S);
S_prew_std = std(epi_S);
E_prew_mean = mean(epi_E);
E_prew_std = std(epi_E);
I_prew_mean = mean(epi_I);
I_prew_std = std(epi_I);
R_prew_mean = mean(epi_R);
R_prew_std = std(epi_R);
    
    
%% Make Figure 4 prewiring vs baseline dynamics
figure(4);
set(gcf, 'Position',  [400, 100, 700,600])% set position, width and height of plot
steps = 1 : (10*24*6) : N_steps;
steps_freq = 1 : (24*6) : N_steps;
color = [0.0,0.4,1.0; 1.0,0.6,0.2; 1.0,0.0,0.2; 0.0,0.8,0.2];


subplot(2,1,1); %baseline

hold all
% % plot S
patch([T_vec(steps) fliplr(T_vec(steps))],...
    [S_bas_mean(steps)+S_bas_std(steps) fliplr(S_bas_mean(steps)-S_bas_std(steps))],...
    color(1,:),'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
plot(T_vec(steps_freq),S_bas_mean(steps_freq),'LineWidth',2,'Color', color(1,:));
% % plot E
patch([T_vec(steps) fliplr(T_vec(steps))],...
    [E_bas_mean(steps)+E_bas_std(steps) fliplr(E_bas_mean(steps)-E_bas_std(steps))],...
    color(2,:),'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
plot(T_vec(steps_freq),E_bas_mean(steps_freq),'LineWidth',2,'Color', color(2,:));

% % plot I
patch([T_vec(steps) fliplr(T_vec(steps))],...
    [I_bas_mean(steps)+I_bas_std(steps) fliplr(I_bas_mean(steps)-I_bas_std(steps))],...
    color(3,:),'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
plot(T_vec(steps_freq),I_bas_mean(steps_freq),'LineWidth',2,'Color', color(3,:));
% % plot R

patch([T_vec(steps) fliplr(T_vec(steps))],...
    [R_bas_mean(steps)+R_bas_std(steps) fliplr(R_bas_mean(steps)-R_bas_std(steps))],...
    color(4,:),'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
plot(T_vec(steps_freq), R_bas_mean(steps_freq),'LineWidth',2,'Color', color(4,:));

hold off
legend('S','E','I','R','Location','west','box','off');
xlabel('Time [days]');
ylabel('Number of nodes');
title('No interventions');
ylim([0,200]);
xlim([0,100]);
set(gca, 'FontSize',16);
xlabel('Time [days]');
box on

    
subplot(2,1,2); %prewire 
    
hold all
%plot S
patch([T_vec(steps) fliplr(T_vec(steps))],...
    [S_prew_mean(steps)+S_prew_std(steps) fliplr(S_prew_mean(steps)-S_prew_std(steps))],...
    color(1,:),'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
plot(T_vec(steps_freq), S_prew_mean(steps_freq),'LineWidth',2,'Color', color(1,:));
%plot E
patch([T_vec(steps) fliplr(T_vec(steps))],...
    [E_prew_mean(steps)+E_prew_std(steps) fliplr(E_prew_mean(steps)-E_prew_std(steps))],...
    color(2,:),'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
plot(T_vec(steps_freq), E_prew_mean(steps_freq),'LineWidth',2,'Color', color(2,:));

% plot I
patch([T_vec(steps) fliplr(T_vec(steps))],...
    [I_prew_mean(steps)+I_prew_std(steps) fliplr(I_prew_mean(steps)-I_prew_std(steps))],...
    color(3,:),'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');

plot(T_vec(steps_freq),I_prew_mean(steps_freq),'LineWidth',2,'Color', color(3,:));

% plot R
patch([T_vec(steps) fliplr(T_vec(steps))],...
    [R_prew_mean(steps)+R_prew_std(steps) fliplr(R_prew_mean(steps)-R_prew_std(steps))],...
    color(4,:),'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');

plot(T_vec(steps_freq), R_prew_mean(steps_freq),'LineWidth',2,'Color', color(4,:));

hold off
legend('S','E','I','R','Location','west','box','off');
title('Prevention');
xlabel('Time [days]');
ylabel('Number of nodes');
ylim([0,200]);
xlim([0,100]);   
set(gca, 'FontSize',16);   
box on
