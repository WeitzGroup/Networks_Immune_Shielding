clc
clear
addpath Data
addpath Functions

%Time parameters
Dt = 1/(24*6); % Time interval in days (10 min)
T = 100; % Total time of simulation (days)
pars.T_vec = 0:Dt:T; % Time vector in unit of days
pars.N_steps = size(pars.T_vec, 2);  % of time steps
% Epidemic parameters
beta = 0.5; % Transmission rate (per day)
gamma_e = 1/2; % Rate from exposed to infected (per day)
gamma_i = 1/6; % Recovery rate (per day)

% Transition probability parameters
pars.P_contact = rate2prob(beta, Dt); % Contact prob. at each time step
pars.P_EI = rate2prob(gamma_e, Dt);
pars.P_IR = rate2prob(gamma_i, Dt);
pars.P_isolation = 1;
n = 200;
pars.n = n;
zero_vector = zeros(n, 1);
pars.zero_vector = zero_vector;

%Id for hcws and patients
is_hcw = zero_vector;
is_hcw(hcws_id) = 1;
is_pat = 1 - is_hcw;
is_hcw = logical(is_hcw);
is_pat = logical(is_pat);
save('Data/hcw_pat_id', 'is_hcw', 'is_pat')
save('Data/epi_params.mat', 'pars')
