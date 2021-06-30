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

% Number of iterations and initial conditions
inf_0 = 1;
R_init = [20, 80]; % initial number of immunized indiv.
n_runs = 500;

% Vectors
epi = zeros(5, N_steps); % 1 - S (0), 2 - E (1),3 - I(2), 4 - R(3), 5 - ISOLATED(4)

%% BASELINE CASE
Rinf_bas = zeros(n_runs, length(R_init));
for i = 1:length(R_init)
    for j = 1:n_runs
    rec_0 = R_init(i);
    node_status = initial_cond(inf_0, rec_0, pars);
    epi = zeros(4, N_steps);
    for i_t = 1:N_steps
        node_status = SEIR_stochastic_fct(adj, node_status, pars); %SEIR dynamics
        epi(:,i_t) = type_to_count(node_status); % Count each type
    end
    Rinf_bas(j, i) = epi(4,end);
    end
end

%% PREWIRING CASE
Rinf_prew = zeros(n_runs, length(R_init));
for i = 1:length(R_init)   
    for j = 1:n_runs
    rec_0 = R_init(i);
    node_status = initial_cond(inf_0, rec_0, pars);
    epi = zeros(4, N_steps);
    %prewire the network
    adj_new = pre_rewire(adj, node_status, is_pat, is_hcw);
    
    for i_t = 1:N_steps
        node_status = SEIR_stochastic_fct(adj_new, node_status, pars); %SEIR dynamics
        epi(:,i_t) = type_to_count(node_status); % Count each type
    end
    Rinf_prew(j, i) = epi(4, end); 
    end   
end
%% ISOLATION CASE
Rinf_iso = zeros(n_runs, length(R_init));

for i = 1:length(R_init)  
    for j = 1:n_runs
    rec_0 = R_init(i);
    node_status = initial_cond(inf_0, rec_0, pars);
    epi = zeros(4, N_steps);    
    
    for i_t = 1:N_steps 
        node_status = SEIR_stochastic_fct(adj, node_status, pars);

        if mod(i_t, isolation_freq) == 0
             node_status = isolate_hcw(node_status, is_hcw, pars); % isolation intervention
        end
        epi(:,i_t) = type_to_count(node_status); % Count each type
    end
    
    % Run de SEIR model over 100 days and isolate infected HCW every week  
    Rinf_iso(j, i) = epi(4, end);

    end  
end
%% REWIRING CASE
Rinf_rew = zeros(n_runs, length(R_init));

for i = 1:length(R_init)
    for j = 1:n_runs
        adj_new = adj;
        rec_0 = R_init(i);
        node_status = initial_cond(inf_0, rec_0, pars);
        epi = zeros(4, N_steps);
        
        for i_t = 1:N_steps    
            node_status = SEIR_stochastic_fct(adj_new, node_status, pars); %SEIR dynamics
            
            if rem(i_t, rewire_freq) == 0
                adj_new = rewire_all(adj_new, node_status, is_pat, is_hcw);
            end
            epi(:, i_t) = type_to_count(node_status);
        end
        Rinf_rew(j, i)=epi(4, end);   
    end
end
%% REWIRING + ISOLATION CASE
Rinf_rewiso = zeros(n_runs, length(R_init));

for i = 1:length(R_init)
    for j = 1:n_runs
        adj_new = adj;
        rec_0 = R_init(i);
        node_status = initial_cond(inf_0, rec_0, pars);
        epi = zeros(4, N_steps);
        
        for i_t = 1:N_steps
            node_status = SEIR_stochastic_fct(adj_new, node_status, pars);
            
            if mod(i_t, isolation_rewiring_freq) == 0
                node_status = isolate_hcw(node_status, is_hcw, pars); % first do isolation of infected hcws
                adj_new = rewire_all(adj_new, node_status, is_pat, is_hcw); % then do rewiring
            end
            epi(:, i_t) = type_to_count(node_status);
        end
        Rinf_rewiso(j, i) = epi(4, end);
    end
end
%% Supplementary Figure 2 - Comparison of probability of outbreak for different interventions and immunity levels
figure(5);
set(gcf, 'Position',  [200, 200, 1200, 450])% set position, width and height of plot
hold all;
colorpalette=[0.7,0.7,0.7;0.9,0.6,0.6;0.3,0.8,0.2;0.1,0.2,0.8;0.6,0.1,0.7];

subplot(1,2,1); %Immunized 10%
%5 stands for no intervention, prew, iso, rew, iso+rew
P_outbreak = zeros(n - R_init(1), 5);
temp_sum = R_init(1) + inf_0;
for x = 0:(n - temp_sum)
    P_outbreak(x+1, 1) = sum((Rinf_bas(:,1) - temp_sum)>=x)./n_runs;
end

for x=0:(n - temp_sum)
    P_outbreak(x+1, 2)=sum((Rinf_prew(:,1) - temp_sum)>=x)./n_runs;
end

for x=0:(n - temp_sum)
    P_outbreak(x+1, 3)=sum((Rinf_iso(:,1) - temp_sum)>=x)./n_runs;
end

for x=0:(n - temp_sum)
    P_outbreak(x+1, 4)=sum((Rinf_rew(:,1) - temp_sum)>=x)/n_runs;
end

for x=0:(n - temp_sum)
    P_outbreak(x+1, 5)=sum((Rinf_rewiso(:,1) - temp_sum)>=x)./n_runs;
end

hold all;
for y = 1:5    
    plot(linspace(0, 1, n - R_init(1)), P_outbreak(:,y), 'LineWidth',2,...
        'Color',colorpalette(y, :));
    xlabel('Total infected [frac.]');
    ylabel('Prob. outbreak')
    legend('no interventions','prewiring','isolation (w)','rewiring (w)',...
        'isolation+rewiring (w)','box','off')
    title(strcat([num2str(100*R_init(1)./n),'% vaccination']));
    set(gca, 'FontSize', 16)     
    box on

end
hold off;

subplot(1,2,2);% Immunized 40%
P_outbreak = zeros(n - R_init(2), 5);
temp_sum = R_init(2) + inf_0;

for x=0:(n-temp_sum)
    P_outbreak(x+1,1)=sum((Rinf_bas(:,2) - temp_sum)>=x)./n_runs;
end
for x=0:(n-temp_sum)
    P_outbreak(x+1,2)=sum((Rinf_prew(:,2) - temp_sum)>=x)./n_runs;
end
for x=0:(n-temp_sum)
    P_outbreak(x+1,3)=sum((Rinf_iso(:,2) - temp_sum)>=x)./n_runs;
end
for x=0:(n-temp_sum)
    P_outbreak(x+1,4)=sum((Rinf_rew(:,2) - temp_sum)>=x)/n_runs;
end
for x=0:(n-temp_sum)
    P_outbreak(x+1,5)=sum((Rinf_rewiso(:,2) - temp_sum)>=x)./n_runs;
end
hold all;

for y=1:5    
    plot(linspace(0, 1, n - R_init(2)), P_outbreak(:,y),'LineWidth',2,...
        'Color', colorpalette(y,:));
   
    xlabel('Total infected [frac.]');
    ylabel('Prob. outbreak');
    legend('no interventions','prewiring','isolation (w)','rewiring (w)',...
        'isolation+rewiring (w)','box','off');
    title(strcat([num2str(100*R_init(2)./n),'% vaccination']));
    set(gca, 'FontSize', 16);     
    box on

end
hold off;
