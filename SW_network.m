function SW_exp = SW_network(varargin)
addpath ../test_rbm_data/
addpath ../power-law-estimator/

% (1) dE/dI must be in (0.086,0.73) when LE/LI = 16 for a bi-stable process

%% Model parameters
nr_events = 1e5;
save_workspace = 0;
save_entire_workspace = 0;


which_process = 'two-state';
% which_process = 'three-state';

% (ER) Erdos-Renyi or (SF) Scale-free or (SF-Chung-Lu) or Star
wake_graph = 'SF-Chung-Lu';
sleep_graph = 'ER';
% Size - parameter
nW = 100; parW = 4.00;
nS = 100; parS = 2.20;

% Probability of forming an inhibitory edge when inh rule is uniform
inhParW = 10;
inhParS = 15;

% Firing rates
lambda_i = 0.001;
lambda_e = 0.013;

% p -value estimation of power-law fits to bouts
fit_PL = 1;
nr_reps = 25;
display_p_val_stuff = 0;
need_apKS = 0;

compute_domains = 1;
compute_bouts = 1;

use_sample_experiment = 0;
% ------------------------------------------------------------------------
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
            inhParW = varargin{i+1}(3);
        case 'S'
            nS = varargin{i+1}(1);
            parS = varargin{i+1}(2);
            inhParS = varargin{i+1}(3);
        case 'save_entire_workspace'
            save_entire_workspace = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
WS_inh_rule = 'ER-like';    % Uniform or Scale-free
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

%% Network configuration generation
use_undirected_inhibitory_edges_accross_graphs = 0;

% Firing and relaxation rates
lambda_b = 0.003;
ri = 0.002;
rE = 0.005;
if strcmp(which_process,'two-state')
    ri = 0;
    rE = 0;
end

% First reaction (original code) or 
% next reaction algorithm (gillespie's algorithm)
use_gillespie_algorithm = 1;

if ~use_sample_experiment
    % adjacency matrices and mean degrees
    % WW(i,j)=1 means, i-th wake node is linked to j-th wake node, 
    % and vice versa
    % dW is the mean number of outgoing excitatory degree of a wake node
    [dW,WW] = generate_random_graph(wake_graph, nW, parW);
    [dS,SS] = generate_random_graph(sleep_graph, nS, parS);
    
    % Adjacency matrix generation for inhibition and mean inh deg
    % WS(i,j)=1 means, i-th wake node inhibits j-th sleep node,
    % SW(i,j)=1 means, i-th sleep node inhibits j-th wake node
    % dIW is the mean number of outgoing inhibitory degree of a wake node
    
    if use_undirected_inhibitory_edges_accross_graphs
        if nW*inhParW == nS*inhParS
            [dIW, WS] = generate_inhibitory_links([nW nS], inhParW, ...
                WW, SS, WS_inh_rule);
            SW = WS';
            dIS = sum(sum(SW))/nS;
        else
            error(['Number of wake-to-sleep inhibitory edges ' ...
                '(%ix%1.2f=%i) must be equal to the number of ' ...
                'sleep-to-wake inhibitory edges (%ix%1.2f=%i)!'], ...
                nW, inhParW, nW*inhParW, nS, inhParS, nS*inhParS);
        end
    else
        [dIW, WS] = generate_inhibitory_links([nW nS], inhParW, ...
            WW, SS, WS_inh_rule);
        [dIS, SW] = generate_inhibitory_links([nS nW], inhParS, ...
            SS, WW, SW_inh_rule);
    end
    
    % Set simulation title
    sim_title = ['W=' wake_graph '(' num2str(nW) ',' ...
        num2str(dW, '%2.2f') ',' num2str(dIW, '%2.2f') ...
        '), S=' sleep_graph '(' num2str(nS) ',' num2str(dS, '%2.2f') ...
        ',' num2str(dIS, '%2.2f') ')' ', # of events=' ...
        num2str(nr_events,'%1.0e, ') which_process];
elseif use_sample_experiment
    %     save('sample_experiments.mat','sample_exps');
    load('sample_experiments.mat','sample_exps');
    nr_events = sample_exps(which_sample_experiment).nr_events;
    sim_title = sample_exps(which_sample_experiment).sim_title;
    WW = sample_exps(which_sample_experiment).W.excitation_matrix;
    dW = sample_exps(which_sample_experiment).W.average_excitation_degree;
    SS = sample_exps(which_sample_experiment).S.excitation_matrix;
    dS = sample_exps(which_sample_experiment).S.average_excitation_degree;
    wake_graph = sample_exps(which_sample_experiment).W.graph;
    sleep_graph = sample_exps(which_sample_experiment).S.graph;
    nW = sample_exps(which_sample_experiment).W.size;
    nS = sample_exps(which_sample_experiment).S.size;
    SW = sample_exps(which_sample_experiment).S.inhibition_matrix;
    dIS = sample_exps(which_sample_experiment).S.average_inhibition_degree;
    WS = sample_exps(which_sample_experiment).W.inhibition_matrix;
    dIW = sample_exps(which_sample_experiment).W.average_inhibition_degree;   
end

%% Model simulation
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
Wnodes_permuted = randperm(nW);
Snodes_permuted = randperm(nS);
W(Wnodes_permuted(1:round(nW*WE_init))) = 1;
W(Wnodes_permuted(end-round(nW*WI_init)+1:end)) = -1;
S(Snodes_permuted(1:round(nS*SE_init))) = 1;
S(Snodes_permuted(end-round(nS*SI_init)+1:end)) = -1;

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
        fprintf(['%%%i of the simulations is completed in ' ...
            '%3.2f minutes\n'], round((100*i)/nr_events), ...
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

%% Mean Field ODE stable equilibrium solutions
[~, ~] = MF_model(WI(1)/nW, WE(1)/nW, SI(1)/nS, SE(1)/nS, ...
    dW, dS, dIW, dIS, lambda_i, lambda_b, lambda_e, ri, rE, ...
    t, WE, SE, [nW nS], plot_simulation, sim_title, which_process);

%% Plot heat map (rows represent WE+1, columns represetnt SE+1)
A = accumarray([(WE+1) (SE+1)], 1, [nW+1 nS+1]);
if plot_heat_map
    figure, surf(A); title(sim_title); axis tight;
    xlabel('# of excited SLEEP nodes');
    ylabel('# of excited WAKE nodes');
end

%% %% Record simulation results in a single structure
SW_exp.sim_title = sim_title;
SW_exp.nr_events = nr_events;

SW_exp.W.graph = wake_graph;
SW_exp.W.size = nW;
SW_exp.W.excitation_matrix = WW;
SW_exp.W.average_excitation_degree = dW;
SW_exp.W.inhibition_matrix = WS;
SW_exp.W.average_inhibition_degree = dIW;

SW_exp.S.graph = sleep_graph;
SW_exp.S.size = nS;
SW_exp.S.excitation_matrix = SS;
SW_exp.S.average_excitation_degree = dS;
SW_exp.S.inhibition_matrix = SW;
SW_exp.S.average_inhibition_degree = dIS;

SW_exp.WE_SE.density = A;

%% Work duration matrix using transition path based state assignment method
% through calculate_ints_tbp
if compute_domains
    which_method = 'new';
    if strcmp(which_method, 'old')
        [duration_matrix, wake_active_domain, sleep_active_domain, ...
            ~, ~] = calculate_ints_tpb(...
            [nW nS], nr_events, sim_title, t, WE, SE, A);
    elseif strcmp(which_method, 'new')
        [wake_active_domain, sleep_active_domain] = ...
            compute_wake_sleep_domains_NEW([nW nS], sim_title, A);
        duration_matrix = compute_duration_matrix(t, WE, SE, ...
            wake_active_domain, sleep_active_domain, nr_events, ...
            sim_title, ...
            'plot_pdfs_for_all_intervals', plot_pdfs_for_all_intervals);
    end
else
    duration_matrix = NaN;
    wake_active_domain = NaN;
    sleep_active_domain = NaN;
end

SW_exp.duration_matrix = duration_matrix;
SW_exp.WE_SE.wake_active_domain = wake_active_domain;
SW_exp.WE_SE.sleep_active_domain = sleep_active_domain;

%% Compute bouts
if compute_bouts
    [sleep_bouts_prev_trans, wake_bouts_prev_trans, ~] = ...
        interval_to_bout(duration_matrix, ...
        'transition intervals are added to previous bout');
    [sleep_intervals, wake_intervals, transition_intervals] = ...
        interval_to_bout(duration_matrix, 'intervals are bouts');
    [sleep_bouts_half_trans, wake_bouts_half_trans, ~] = ...
        interval_to_bout(duration_matrix, ...
        ['transition intervals are cut in half and ' ...
        'added to closest wake/sleep interval']);
else
    sleep_bouts_prev_trans = [];
    wake_bouts_prev_trans = [];
    sleep_intervals = [];
    wake_intervals = [];
    transition_intervals = [];
    sleep_bouts_half_trans = [];
    wake_bouts_half_trans = [];
end

SW_exp.intervals.wake = wake_intervals;
SW_exp.intervals.sleep = sleep_intervals;
SW_exp.intervals.transition = transition_intervals;

SW_exp.bouts.wake.no_trans.data = wake_intervals;
SW_exp.bouts.wake.prev_trans.data = wake_bouts_prev_trans;
SW_exp.bouts.wake.half_trans.data = wake_bouts_half_trans;
SW_exp.bouts.sleep.no_trans.data = sleep_intervals;
SW_exp.bouts.sleep.prev_trans.data = sleep_bouts_prev_trans;
SW_exp.bouts.sleep.half_trans.data = sleep_bouts_half_trans;

%% Fit power-law
if fit_PL
    PL_sbpt = find_pl_fit_with_adapt_slope_detect(...
        sleep_bouts_prev_trans, ...
        'data_title', 'Sleep bouts (prev-trans)', ...
        'display_p_val_stuff', display_p_val_stuff, ...
        'nr_reps', nr_reps, ...
        'plot_best_pl_fit', 0, 'need_apKS', need_apKS);
    PL_wbpt = find_pl_fit_with_adapt_slope_detect(...
        wake_bouts_prev_trans, ...
        'data_title', 'Wake bouts (prev-trans)', ...
        'display_p_val_stuff', display_p_val_stuff, ...
        'nr_reps', nr_reps, ...
        'plot_best_pl_fit', 0, 'need_apKS', need_apKS);
    PL_sbnt = find_pl_fit_with_adapt_slope_detect(...
        sleep_intervals, ...
        'data_title', 'Sleep bouts (no-trans)', ...
        'display_p_val_stuff', display_p_val_stuff, ...
        'nr_reps', nr_reps, ...
        'plot_best_pl_fit', 0, 'need_apKS', need_apKS);
    PL_wbnt = find_pl_fit_with_adapt_slope_detect(...
        wake_intervals, ...
        'data_title', 'Wake bouts (no-trans)', ...
        'display_p_val_stuff', display_p_val_stuff, ...
        'nr_reps', nr_reps, ...
        'plot_best_pl_fit', 0, 'need_apKS', need_apKS);
    PL_sbht = find_pl_fit_with_adapt_slope_detect(...
        sleep_bouts_half_trans, ...
        'data_title', 'Sleep bouts (half-trans)', ...
        'display_p_val_stuff', display_p_val_stuff, ...
        'nr_reps', nr_reps, ...
        'plot_best_pl_fit', 0, 'need_apKS', need_apKS);
    PL_wbht = find_pl_fit_with_adapt_slope_detect(...
        wake_bouts_half_trans, ...
        'data_title', 'Wake bouts (half-trans)', ...
        'display_p_val_stuff', display_p_val_stuff, ...
        'nr_reps', nr_reps, ...
        'plot_best_pl_fit', 0, 'need_apKS', need_apKS);
else
    PL_sbpt = [];
    PL_wbpt = [];
    PL_sbnt = [];
    PL_wbnt = [];
    PL_sbht = [];
    PL_wbht = [];
end

SW_exp.bouts.wake.no_trans.PL = PL_wbnt;
SW_exp.bouts.wake.prev_trans.PL = PL_wbpt;
SW_exp.bouts.wake.half_trans.PL = PL_wbht;
SW_exp.bouts.sleep.no_trans.PL = PL_sbnt;
SW_exp.bouts.sleep.prev_trans.PL = PL_sbpt;
SW_exp.bouts.sleep.half_trans.PL = PL_sbht;

%% Fit exponenential
if fit_exp_to_tail
    exp_tail_sbpt = fit_exponential_to_tail(sleep_bouts_prev_trans, ...
        'plot_stuff', 0);
    exp_tail_wbpt = fit_exponential_to_tail(wake_bouts_prev_trans, ...
        'plot_stuff', 0);
    exp_tail_sbnt = fit_exponential_to_tail(sleep_intervals, ...
        'plot_stuff', 0);
    exp_tail_wbnt = fit_exponential_to_tail(wake_intervals, ...
        'plot_stuff', 0);
    exp_tail_sbht = fit_exponential_to_tail(sleep_bouts_half_trans, ...
        'plot_stuff', 0);
    exp_tail_wbht = fit_exponential_to_tail(wake_bouts_half_trans, ...
        'plot_stuff', 0);
else
    exp_tail_sbpt = [];
    exp_tail_wbpt = [];
    exp_tail_sbnt = [];
    exp_tail_wbnt = [];
    exp_tail_sbht = [];
    exp_tail_wbht = [];
end

SW_exp.bouts.wake.prev_trans.exp_tail = exp_tail_wbpt;
SW_exp.bouts.wake.no_trans.exp_tail = exp_tail_wbnt;
SW_exp.bouts.wake.half_trans.exp_tail = exp_tail_wbht;
SW_exp.bouts.sleep.no_trans.exp_tail = exp_tail_sbnt;
SW_exp.bouts.sleep.prev_trans.exp_tail = exp_tail_sbpt;
SW_exp.bouts.sleep.half_trans.exp_tail = exp_tail_sbht;

%% Calculate drift-diffusion coefficients
if calculate_drift_diffusion
    [F,D] = calc_drift_diffusion_potential('2D', [nW nS], ...
        sim_title, t, WE, SE);
    [F_WE,D_WE] = calc_drift_diffusion_potential('1D', nW, ...
        sim_title, t, WE, 'WE', 'plot_DDE_stuff', 0);
    [F_SE,D_SE] = calc_drift_diffusion_potential('1D', nS, ...
        sim_title, t, SE, 'SE', 'plot_DDE_stuff', 0);
else
    F = zeros(nW+1,nS+1,2,1);
    D = zeros(2,2,nW+1,nS+1);
    F_WE = zeros(nW+1,1);    
    D_WE = zeros(nW+1,1);
    F_SE = zeros(nS+1,1);    
    D_SE = zeros(nS+1,1);
end

SW_exp.DDE.F = F;
SW_exp.DDE.D = D;
SW_exp.DDE.F_WE = F_WE;
SW_exp.DDE.D_WE = D_WE;
SW_exp.DDE.F_SE = F_SE;
SW_exp.DDE.D_SE = D_SE;

SW_exp.use_gillespie_algorithm = use_gillespie_algorithm;

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

% Plot PL
if fit_PL
    plot_sleep_wake_bouts_PL_results(SW_exp);
end

% Plot bout distributions
if plot_bout_dist
    N_sleep = hist(sleep_bouts, ...
        round((max(sleep_bouts)-min(sleep_bouts))/10));
    N_wake = hist(wake_bouts, ...
        round((max(wake_bouts)-min(wake_bouts))/10));
    
    figure,
    subplot(2,2,1);
    loglog(N_sleep,'b-');
    title('straight ==> Sleep bouts are power-law');
    subplot(2,2,2);
    semilogy(N_sleep,'b-');
    title('straight ==> Sleep bouts are exponential');
    subplot(2,2,3);
    loglog(N_wake,'b-');
    title('straight ==> Wake bouts are power-law');
    subplot(2,2,4);
    semilogy(N_wake,'b-');
    title('straight ==> Wake bouts are exponential');
end

% Plot degree distributions
if plot_degree_dist
    figure,
    subplot(2,2,1),
    hist(sum(WW),10);
    title('Wake cluster degree distribution');
    subplot(2,2,2),
    hist(sum(SS),10);
    title('Sleep cluster degree distribution');
    subplot(2,2,3),
    hist(sum(WS,2),10);
    title('Wake to Sleep inhibitory degree distribution');
    subplot(2,2,4),
    hist(sum(SW,2),10);
    title('Sleep to Wake inhibitory degree distribution');
end

% Give warning of possible non-oscillation stochastic process when
% transition intervals take up too much time
% total_trans_time = sum(duration_matrix(duration_matrix(:,2)==0,1));
% if total_trans_time/t(end) > 0.10
%     warning('Fatih: Total transition time if too much. Stochastic 
% process may not be oscillating.');

fprintf('\nENTIRE SIMULATION TOOK %3.2f minutes\n',toc(tSIM)/60);

end