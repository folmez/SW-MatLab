function SW_exp = SW_network(varargin)
%% Model parameters
nr_events = 1e6;
save_workspace = 0;
save_entire_workspace = 0;

which_process = 'two-state';    % 'two-state' or 'three-state'

% (ER) Erdos-Renyi or (SF) Scale-free or (SF-Chung-Lu) or Star
wake_graph = 'ER';
sleep_graph = 'ER';
% Size - parameter
nW = 100; parW = 3.00;
nS = 100; parS = 3.00;
% Average inhibitory degree when inhibitin rule is ER-like
inh_parW = 8.;
inh_parS = 8.;
% Note: dE/dI must be in (0.086,0.73) when LE/LI = 16 for a bi-stable
% process

% Firing rates
lambda_i = 0.001;
lambda_e = 0.016;

% p -value estimation of power-law fits to bouts
fit_PL = 1;
nr_reps = 25;
display_p_val_stuff = 0;
need_apKS = 0;

compute_domains = 1;
compute_bouts = 1;

use_sample_experiment = 0;
which_sample_experiment = 0;

i=1;
while i<=length(varargin),
    switch varargin{i},
        case 'use_sample_experiment'
            use_sample_experiment = varargin{i+1}(1);
            which_sample_experiment = varargin{i+1}(2);
        case 'save_workspace',      save_workspace = varargin{i+1};
        case 'which_process',       which_process = varargin{i+1};
        case 'nr_events',           nr_events = varargin{i+1};
        case 'compute_domains',     compute_domains = varargin{i+1};
        case 'compute_bouts',       compute_bouts = varargin{i+1};
        case 'firing_rates'
            lambda_i = varargin{i+1}(1);
            lambda_e = varargin{i+1}(2);
        case 'fit_PL'              
            fit_PL = varargin{i+1}(1);
            nr_reps = varargin{i+1}(2);
            display_p_val_stuff = varargin{i+1}(3);
            need_apKS = varargin{i+1}(4);
        case 'wake_graph',          wake_graph = varargin{i+1};
        case 'sleep_graph',         sleep_graph = varargin{i+1};
        case 'W'
            nW = varargin{i+1}(1);
            parW = varargin{i+1}(2);
            inh_parW = varargin{i+1}(3);
        case 'S'
            nS = varargin{i+1}(1);
            parS = varargin{i+1}(2);
            inh_parS = varargin{i+1}(3);
        case 'save_entire_workspace'
            save_entire_workspace = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

WS_inh_rule = 'ER-like';
SW_inh_rule = 'ER-like';

display_event_sim_summary = 0;

plot_simulation = 1;
plot_bout_dist = 0;
plot_heat_map = 1;
plot_degree_dist = 0;
plot_pdfs_for_all_intervals = 1;
plot_sample_trajectories = 0;

calculate_drift_diffusion = 0;

fit_exp_to_tail = 1;

% Firing and relaxation rates
lambda_b = 0.003;
if strcmp(which_process,'two-state')
    ri = 0;
    rE = 0;
elseif strcmp(which_process,'three-state')
    ri = 0.002;
    rE = 0.005;
end

% First reaction (original code) or 
% next reaction algorithm (gillespie's algorithm)
use_gillespie_algorithm = 1;

use_undirected_inhibitory_edges_accross_graphs = 0;

%% Network configuration generation
[nr_events, sim_title, WW, dW, SS, dS, wake_graph, sleep_graph, ...
    nW, nS, SW, dIS, WS, dIW] = ...
    generate_network_configuration(...
    use_sample_experiment, which_sample_experiment, ...
    wake_graph, [nW, parW, inh_parW], ...
    sleep_graph, [nS, parS, inh_parS], ...
    WS_inh_rule, SW_inh_rule, ...
    use_undirected_inhibitory_edges_accross_graphs, ...
    nr_events, which_process, plot_degree_dist);

%% Model simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial number of excited and inhibited nodes in W and S
WE_init = rand;
SE_init = rand;
if strcmp(which_process,'three-state')
    WI_init = rand*(1-WE_init);
    SI_init = rand*(1-SE_init);
elseif strcmp(which_process,'two-state')
    WI_init = 1-WE_init;
    SI_init = 1-SE_init;
end

% Initial states of nodes: -1 for I, 0 for B, 1 for E
W = zeros(nW,1); S = zeros(nS,1);
W_nodes_permuted = randperm(nW);
S_nodes_permuted = randperm(nS);
W(W_nodes_permuted(1:round(nW*WE_init))) = 1;
W(W_nodes_permuted(end-round(nW*WI_init)+1:end)) = -1;
S(S_nodes_permuted(1:round(nS*SE_init))) = 1;
S(S_nodes_permuted(end-round(nS*SI_init)+1:end)) = -1;

% Sizes of compartmens at every time
WI = zeros(nr_events,1); WE = zeros(nr_events,1);
SI = zeros(nr_events,1); SE = zeros(nr_events,1);
t = zeros(nr_events,1);

% Initiation
WI(1) = sum(W==-1); WE(1) = sum(W==1);
SI(1) = sum(S==-1); SE(1) = sum(S==1);

tSIM = tic;
tDisplay = tic;
for i=2:nr_events
    [W, S, wait_time] = update_SW_network(i, W, S, WI, WE, ...
        SI, SE, WW, SS, WS, SW, [nW nS], lambda_i, lambda_b, lambda_e, ...
        ri, rE, use_gillespie_algorithm, display_event_sim_summary, ...
        which_process);

    % Update variables and time
    WI(i) = sum(W==-1);
    WE(i) = sum(W==1);
    SI(i) = sum(S==-1);
    SE(i) = sum(S==1);
    t(i) = t(i-1)+wait_time;

    % Progress report
    if toc(tDisplay)>6
        fprintf(['%%%2.2f of the simulations is completed in ' ...
            '%3.2f minutes\n'], (100*i)/nr_events, ...
            toc(tSIM)/60);
        tDisplay = tic;
    end    
    
    % Show surface at an intermediate time to confirm bistability for
    % single runs
    if i == 1e5
        A = accumarray([(WE(1:min(nr_events,1e5))+1) ...
            (SE(1:min(nr_events,1e5))+1)], 1, [nW+1 nS+1]);
        figure
        surf(A);
        axis tight;
        title('Press a key to continue if bi-stable...', 'FontSize', 20);
        drawnow;
        pause;
    end
end
fprintf('Competitive graph model simulation took %3.2f minutes\n', ...
    toc(tSIM)/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mean Field ODE stable equilibrium solutions
[~, ~] = MF_model(WI(1)/nW, WE(1)/nW, SI(1)/nS, SE(1)/nS, ...
    dW, dS, dIW, dIS, lambda_i, lambda_b, lambda_e, ri, rE, ...
    t, WE, SE, [nW nS], plot_simulation, sim_title, which_process);

%% Plot heat map (rows represent WE+1, columns represetnt SE+1)
A = accumarray([(WE+1) (SE+1)], 1, [nW+1 nS+1]);
if plot_heat_map
    figure,
    surf(A);
    title(sim_title);
    axis tight;
    xlabel('# of excited SLEEP nodes');
    ylabel('# of excited WAKE nodes');
end

%% Compute activity domains using transition path based assignment
[duration_matrix, wake_active_domain, sleep_active_domain] = ...
    compute_activity_domains([nW nS], nr_events, sim_title, t, WE, SE, ...
    A, plot_pdfs_for_all_intervals, compute_domains);

%% Compute bouts
[sleep_bouts_prev_trans, wake_bouts_prev_trans, ...
    sleep_intervals, wake_intervals, transition_intervals, ...
    sleep_bouts_half_trans, wake_bouts_half_trans] = ...
    compute_bout_durations(duration_matrix, compute_bouts);

%% Fit bounded power-law and exponential to SW bouts
[PL_sbpt, PL_wbpt, PL_sbnt, PL_wbnt, PL_sbht, PL_wbht, ...
    exp_tail_sbpt, exp_tail_wbpt, exp_tail_sbnt, ...
    exp_tail_wbnt, exp_tail_sbht, exp_tail_wbht] = ...
    fit_PL_and_EXP_to_SW_bouts(...
    sleep_bouts_prev_trans, wake_bouts_prev_trans, ...
    sleep_intervals, wake_intervals, ...
    sleep_bouts_half_trans, wake_bouts_half_trans, ...
    display_p_val_stuff, nr_reps, fit_PL, fit_exp_to_tail);

%% Calculate drift-diffusion coefficients
[F, D, F_WE, D_WE, F_SE, D_SE] = do_drift_diffusion_analysis(...
    [nW nS], sim_title, t, WE, SE);

%% %% Record simulation results in a single structure
SW_exp = record_sim(sim_title, nr_events, ...
    wake_graph,  nW, WW, dW, WS, dIW, ...
    sleep_graph, nS, SS, dS, SW, dIS, A, ...
    duration_matrix, wake_active_domain, sleep_active_domain, ...
    wake_intervals, sleep_intervals, transition_intervals, ...
    wake_bouts_prev_trans, wake_bouts_half_trans, ...
    sleep_bouts_prev_trans, sleep_bouts_half_trans, ...
    PL_wbnt, PL_wbpt, PL_wbht, PL_sbnt, PL_sbpt, PL_sbht, ...
    exp_tail_wbpt, exp_tail_wbnt, exp_tail_wbht, ...
    exp_tail_sbnt, exp_tail_sbpt, exp_tail_sbht, ...
    F, D, F_WE, D_WE, F_SE, D_SE, ...
    use_gillespie_algorithm);

%% Save workspace
filename = ['saved_workspaces/network_' which_process '_' ...
    wake_graph '_' num2str(nW) '_' num2str(dW) '_' num2str(dIW) '_' ...
    sleep_graph '_' num2str(nS) '_' num2str(dS) '_' num2str(dIS) '_' ...
    datestr(now,'mmmdd_HHMMSS') '.mat'];
if save_entire_workspace
    save([filename(1:end-4) '_ENTIRE.mat']);
elseif save_workspace
    save(filename,'experiment');
end

%% Check outputs
if strcmp(which_process,'two-state') && ...
        ~(prod(WE+WI==nW)==1 && prod(SE+SI==nS)==1)
    error('Sum of excited and inhibited nodes don''t make n!!!')    
end

%% Plot stuff

% Plot sample trajectories
if plot_sample_trajectories
    make_movie_from_bout_trajs = 0;
    what_type_of_bouts = [-1];
    
    min_bout_size = 0;
    max_bout_size = inf;
    nr_bout_trajs = 50;
    pause_after_figure = 1;
    
    if make_movie_from_bout_trajs, close all, end
    
    for i=1:nr_bout_trajs
        bout_index_set = find(duration_matrix(:,1) > min_bout_size & ...
            duration_matrix(:,1) < max_bout_size);
        bout_index = bout_index_set(randperm(length(bout_index_set),1));
        while ~ismember(duration_matrix(bout_index,2),what_type_of_bouts)
            bout_index = bout_index_set( ...
                randperm(length(bout_index_set),1));
        end
        
        figure, hold on;
        plot(wake_active_domain(:,2), wake_active_domain(:,1), 'bo');
        plot(sleep_active_domain(:,2), sleep_active_domain(:,1), 'gd');
        
        for j = 0:2
            bout_begin_time = sum(duration_matrix(1:bout_index+j-1,1));
            bout_end_time = bout_begin_time + ...
                duration_matrix(bout_index+j,1);
            bout_type = duration_matrix(bout_index+j,2);
            
            traj_index_set_min = find(t >= bout_begin_time,1);
            traj_index_set_max = find(t > bout_end_time,1);
            traj_index_set = (traj_index_set_min:traj_index_set_max-1)';
            if bout_type==-1
                str_bout_type = 'Sleep';
                plot(SE(traj_index_set(1)), WE(traj_index_set(1)), 'g*');
                plot(SE(traj_index_set), WE(traj_index_set), 'g.-');
            elseif bout_type==1
                str_bout_type = 'Wake';
                plot(SE(traj_index_set(1)), WE(traj_index_set(1)), 'b*');
                plot(SE(traj_index_set), WE(traj_index_set), 'b.-');
            elseif bout_type==0
                str_bout_type = 'Transition';
                plot(SE(traj_index_set(1)), WE(traj_index_set(1)), 'k*');
                plot(SE(traj_index_set), WE(traj_index_set), 'k.-');
            end
            %plot(SE(traj_index_set(1)), WE(traj_index_set(1)), 'r*');
            %plot(SE(traj_index_set(end)), WE(traj_index_set(end)), 'k*');
        end
        legend(sim_title,'Location','Best');
        axis([0 nS+1 0 nW+1]); xlabel('SE'); ylabel('WE');
        title(['Sample ' str_bout_type ' bout, size = ' ...
            num2str(duration_matrix(bout_index,1))]);
        
        if pause_after_figure, pause; end
    end
end

end