function MFMC_exp = two_ER_sim_MFMC(varargin)
addpath ../test_rbm_data/
addpath ../power-law-estimator/

nr_events = 1e4;
save_workspace = 0;
% ------------------------------------------------------------------------
i=1;
while i<=length(varargin),
    switch varargin{i},
        case 'save_workspace',      save_workspace = varargin{i+1};
        case 'nr_events',           nr_events = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
WE_init = 0.5; WI_init = 0.1;
SE_init = 0.1; SI_init = 0.8;

wake_cluster = 'ER';
sleep_cluster = 'ER';

plot_bout_dist = 1;
plot_degree_dist = 1;
plot_simulation = 1;
plot_heat_map = 1;
calculate_2D_drift_diffusion = 0;
use_gillespie_algorithm = 1; % always 1 
display_event_sim_summary = 1;

N = 100;
mean_deg_tol = 0.005;

% probability of forming an excitatory edge or power of power law distribution for degrees
parW = 0.076;
parS = 0.076; % For Scale-free, mean degree ~ 2*par!!!

WS_inh_rule = 'Uniform';
SW_inh_rule = 'Uniform';

% Probability of forming an inhibitory edge when inh rule is uniform
inhParW = 0.068;
inhParS = 0.068;

ri = 0.002;
rE = 0.005;

lambda_i = 0.001;
lambda_b = 0.003;
lambda_e = 0.016;

% adjacency matrices and mean degrees
% WW(i,j)=1 means, i-th wake node is linked to j-th wake node, and vice versa
% dW is the mean number of outgoing excitatory degree of a wake node
[dW,WW] = generate_random_graph(wake_cluster, N, parW, mean_deg_tol);
[dS,SS] = generate_random_graph(sleep_cluster, N, parS, mean_deg_tol);

% adjacency matrix generation for inhibition and mean inh deg
% WS(i,j)=1 means, i-th wake node inhibits j-th sleep node, SW(i,j)=1 means, i-th sleep node inhibits j-th wake node
% dIW is the mean number of outgoing inhibitory degree of a wake node
[dIW,WS] = generate_inhibitory_links(N, inhParW, WW, SS, WS_inh_rule, mean_deg_tol);
[dIS,SW] = generate_inhibitory_links(N, inhParS, SS, WW, SW_inh_rule, mean_deg_tol);

EdegW = sum(WW)'; % excitatory degree dist of wake cluster
EdegS = sum(SS)'; % excitatory degree dist of sleep cluster
IdegW = sum(WS)'; % inhibitory degree dist of wake cluster
IdegS = sum(WS,2); % inhibitory degree dist of sleep cluster

WE = zeros(nr_events,1); WI = zeros(nr_events,1);
SE = zeros(nr_events,1); SI = zeros(nr_events,1);

t = zeros(nr_events,1); % wait times

% Initiation
WE(1) = WE_init*N; WI(1) = WI_init*N;
SE(1) = SE_init*N; SI(1) = SI_init*N;

sim_title = ['[MFMC] W=' wake_cluster '(' num2str(N) ',' num2str(dW) ...
    ',' num2str(dIW) '), S=' sleep_cluster '(' num2str(N) ',' ...
    num2str(dS) ',' num2str(dIS) ')' ', # of events=' num2str(nr_events)];
%------------------------------------------------------------------
tSIM = tic;
progress_counter = nr_events/100;
for i=2:nr_events
    [WE(i), WI(i), SE(i), SI(i), t(i)] = ...
        update_MFMC_network_process(i, WE(i-1), WI(i-1), ...
        SE(i-1), SI(i-1), t(i-1), EdegW, EdegS, IdegW, IdegS, N, ...
        lambda_i, lambda_b, lambda_e, ri, rE, display_event_sim_summary);
    
    % Progress report
    if i == progress_counter
        fprintf('%%%i of the simulations is completed in %3.2f minutes\n',...
            (100*i)/nr_events, toc(tSIM)/60);
        progress_counter = progress_counter + nr_events/100;
    end
end
fprintf('\nMean-field Markov chain simulation took %3.2f minutes\n', ...
    toc(tSIM)/60);
%------------------------------------------------------------------
% Mean Field ODE stable equilibrium solutions
[t_out, y_out]= MF_model(WI(1)/N, WE(1)/N, SI(1)/N, SE(1)/N, ...
    dW, dS, dIW, dIS, lambda_i, lambda_b, lambda_e, ri, rE, ...
    t, WE, SE, N, plot_simulation, sim_title);
%------------------------------------------------------------------
% Plot heat map (rows represent WE+1, columns represetnt SE+1)
A = accumarray([(WE+1) (SE+1)], 1, [N+1 N+1]);
if plot_heat_map
    figure, surf(A); title(sim_title);
    xlabel('# of excited SLEEP nodes');
    ylabel('# of excited WAKE nodes');
end
%------------------------------------------------------------------
% Work duration matrix using transition path based state assignment method
% through calculate_ints_tbp
[duration_matrix, wake_active_domain, sleep_active_domain, ...
    wake_percentage, sleep_percentage] = calculate_ints_tpb(...
    N, nr_events, sim_title, t, WE, SE, A);
%------------------------------------------------------------------
[wake_bouts, sleep_bouts, ~] = interval_to_bout(...
    duration_matrix,'transition intervals are added to previous bout');
%------------------------------------------------------------------
if calculate_2D_drift_diffusion
    [F,D] = calculate2D_drift_diffusion_potential(N,sim_title,t,WE,SE);
else
    F = zeros(N+1,N+1,2,1);
    D = zeros(2,2,N+1,N+1);
end
%------------------------------------------------------------------
% plot bout distributions
if plot_bout_dist
    N_sleep = hist(sleep_bouts,round((max(sleep_bouts)-min(sleep_bouts))/10));
    N_wake = hist(wake_bouts,round((max(wake_bouts)-min(wake_bouts))/10));
    
    figure,
    subplot(2,2,1); loglog(N_sleep,'b-'), title('straight ==> Sleep bouts are power-law');
    subplot(2,2,2); semilogy(N_sleep,'b-'), title('straight ==> Sleep bouts are exponential');
    subplot(2,2,3); loglog(N_wake,'b-'), title('straight ==> Wake bouts are power-law');
    subplot(2,2,4); semilogy(N_wake,'b-'), title('straight ==> Wake bouts are exponential');
end
%------------------------------------------------------------------
if plot_degree_dist
    figure,
        subplot(2,2,1); hist(sum(WW),10); title('Wake cluster degree distribution');
        subplot(2,2,2); hist(sum(SS),10); title('Sleep cluster degree distribution');
        subplot(2,2,3); hist(sum(WS,2),10); title('Wake to Sleep inhibitory degree distribution');
        subplot(2,2,4); hist(sum(SW,2),10); title('Sleep to Wake inhibitory degree distribution');
end
%------------------------------------------------------------------
% CREATE SINGLE STRUCTURED VARIABLE THAT HOLDS SIMULATION RESULTS
MFMC_exp.sim_title = sim_title;
MFMC_exp.N = N;
MFMC_exp.nr_events = nr_events;
MFMC_exp.wake_cluster = wake_cluster;
MFMC_exp.sleep_cluster = sleep_cluster;
MFMC_exp.WW = WW;
MFMC_exp.SS = SS;
MFMC_exp.WS = WS;
MFMC_exp.SW = SW;
MFMC_exp.dS = dS;
MFMC_exp.dW = dW;
MFMC_exp.dIS = dIS;
MFMC_exp.dIW = dIW;

MFMC_exp.A = A;

MFMC_exp.wake_active_domain = wake_active_domain;
MFMC_exp.sleep_active_domain = sleep_active_domain;
MFMC_exp.wake_percentage = wake_percentage;
MFMC_exp.sleep_percentage = sleep_percentage;

MFMC_exp.duration_matrix = duration_matrix;
MFMC_exp.F = F;
MFMC_exp.D = D;

MFMC_exp.wake_bouts = wake_bouts;
MFMC_exp.sleep_bouts = sleep_bouts;

MFMC_exp.use_gillespie_algorithm = use_gillespie_algorithm;
fprintf('\nENTIRE SIMULATION TOOK %3.2f minutes\n',toc(tSIM)/60);
%------------------------------------------------------------------
if save_workspace
    filename = ['two_cluster_MFMC_' ...
        wake_cluster '_' num2str(dW) '_' num2str(dIW) '_'...
        sleep_cluster '_' num2str(dS) '_' num2str(dIS) '_'...
        datestr(now,'mmmdd_HHMMSS') '.mat'];    
    save(filename,'MFMC_exp');
end
%------------------------------------------------------------------
end