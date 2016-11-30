function [exps, avg_exp] = run_multiple_network_sims(varargin)
addpath ../test_rbm_data/
addpath ../power-law-estimator/

save_workspace = 0;
load_from_file = 0;

nr_events = 1e5;
which_process = 'two-state';
wake_graph = 'ER'; %'SF';   % (ER) Erdos-Renyi or (SF) Scale-free or Star
sleep_graph = 'ER'; %'SF';
nW = 100;
nS = 100;
% probability of forming an excitatory edge or power of power law distribution for degrees
parW = 30;%0.03; %3;
parS = 30;%0.03; %3; % For Scale-free, mean degree ~ 2*par!!!
% Probability of forming an inhibitory edge when inh rule is uniform
inhParW = 8; % 0.045; % doesn't matter when inh rule is not uniform!!!           
inhParS = 8; % 0.045;

% Firing rates
lambda_i = 0.001;
lambda_e = 0.020;

fit_PL = 0;
nr_reps = 25;
display_p_val_stuff = 0;
need_only_KS_method = 0;

plot_avg_density = 1;
% display_PL_results = 0;
compute_domains = 1;
compute_bouts = 1;

use_sample_experiment = 0;
which_sample_experiment = 0;
% ------------------------------------------------------------------------
i=1;
while i<=length(varargin),
    switch varargin{i},
        case 'use_sample_experiment'
            use_sample_experiment = varargin{i+1}(1);
            which_sample_experiment = varargin{i+1}(2);
        case 'load from file',      
            load_from_file = 1;
            filenames = varargin{i+1}';
            [nr_filenames, ~] = size(filenames);
            if nr_filenames==1
                load(varargin{i+1});
            end
            if exist('experiment','var'), exps = experiment; end
        case 'which_process',       which_process = varargin{i+1};
        case 'firing_rates',        
            lambda_i = varargin{i+1}(1);
            lambda_e = varargin{i+1}(2);
        case 'nr_sims',             nr_sims = varargin{i+1};
        case 'compute_domains',     compute_domains = varargin{i+1};
        case 'compute_bouts',       compute_bouts = varargin{i+1};
        case 'save_workspace',      save_workspace = varargin{i+1};
        case 'wake_graph',          wake_graph = varargin{i+1};
        case 'sleep_graph',         sleep_graph = varargin{i+1};
        case 'nr_events',           nr_events = varargin{i+1};
        case 'display_PL_results',  display_PL_results = varargin{i+1};
        case 'fit_PL',              fit_PL = varargin{i+1}(1);
                                    nr_reps = varargin{i+1}(2);
                                    display_p_val_stuff = varargin{i+1}(3);
                                    need_only_KS_method = varargin{i+1}(4);
        case 'plot_avg_density',    plot_avg_density = varargin{i+1};
        case 'p',                   parW = varargin{i+1};
                                    parS = varargin{i+1};
        case 'q',                   inhParW = varargin{i+1};
                                    inhParS = varargin{i+1};
        case 'W'
            nW = varargin{i+1}(1);
            parW = varargin{i+1}(2);
            inhParW = varargin{i+1}(3);
        case 'S'
            nS = varargin{i+1}(1);
            parS = varargin{i+1}(2);
            inhParS = varargin{i+1}(3);
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
plot_individial_densities = 0;

display_mean_PL_results = 1;
display_individual_PL_results = 1;
display_individual_exp_tail_results = 1;
plot_individual_PL_results = 0;

display_progress = 0;
plot_all_bouts = 1;
% ------------------------------------------------------------------------
if ~load_from_file
    % Run multiple competitive graph model simulations
    for i=1:nr_sims
        fprintf('SIMULATION %i/%i\n', i, nr_sims);
        exps(nr_sims+1-i) = two_ER_simulation('save_workspace',0, ...
            'W', [nW parW inhParW], 'wake_graph', wake_graph, ... 
            'S', [nS parS inhParS], 'sleep_graph', sleep_graph, ... 
            'firing_rates', [lambda_i lambda_e], ...
            'nr_events', nr_events, 'which_process', which_process, ...
            'display_progress', display_progress, ...
            'compute_bouts', compute_bouts, ...
            'compute_domains', compute_domains, ...
            'fit_PL', [fit_PL nr_reps display_p_val_stuff need_only_KS_method], ...
            'use_sample_experiment', [use_sample_experiment which_sample_experiment]);
    end
    % Save simulations
    if save_workspace
        dW = exps(1).W.average_excitation_degree;
        dS = exps(1).S.average_excitation_degree;
        dIW = exps(1).W.average_inhibition_degree;
        dIS = exps(1).S.average_inhibition_degree;
        if ~use_sample_experiment
            filename = ['saved_workspaces/mult_network_' which_process '_' ...
                num2str(nr_sims) '_sims_' ...
                wake_graph '_' num2str(nW) '_' num2str(dW) '_' num2str(dIW) '_' ...
                sleep_graph '_' num2str(nS) '_' num2str(dS) '_' num2str(dIS) '_' ...
                'lambda_ratio_' num2str(lambda_e/lambda_i) '_', ...
                datestr(now,'mmmdd_HHMMSS') '.mat'];
        elseif use_sample_experiment
            filename = ['saved_workspaces/mult_network_sample_' ...
                num2str(which_sample_experiment) '_' ...
                num2str(nr_sims) '_sims_' ...
                'lambda_ratio_' num2str(lambda_e/lambda_i) '_', ...
                datestr(now,'mmmdd_HHMMSS') '.mat'];
        end
        save(filename,'exps');
    end
end
% ------------------------------------------------------------------------
if nr_filenames>1
    avg_exp{nr_filenames} = [];
    for k = 1:nr_filenames
        load(filenames{k}, 'exps');
        avg_exp{k} = display_sleep_wake_bouts_PL_results(exps, ...
            1, avg_exp{k});
        avg_exp{k} = display_sleep_wake_bouts_exponential_tail_results(exps, ...
            1, avg_exp{k});
        avg_exp{k}.sim_title = exps(1).sim_title;
        avg_exp{k}.graph_titles.wake = avg_exp{k}.sim_title( ...
            3:find(avg_exp{k}.sim_title==' ')-2);
        avg_exp{k}.graph_titles.sleep = avg_exp{k}.sim_title( ...
            find(avg_exp{k}.sim_title==' ')+3:find(avg_exp{k}.sim_title=='#')-3);
        if plot_all_bouts
            avg_exp{k} = collect_all_bouts_in_avg_exp(exps, avg_exp{k});
            plot_sleep_wake_bouts_PL_results(avg_exp{k}, 1);
        end
    end
    
    display_multiple_SW_coupling_results(avg_exp);    
    display_multiple_SYMMETRIC_SW_coupling_results(avg_exp(1:10));
else
    % Find average experiment
    avg_A = zeros(nW+1,nS+1);
    for i=1:length(exps)
        % All densities
        avg_A = avg_A + exps(i).WE_SE.density;
        if plot_individial_densities
            plot_activity_domains_on_surface_and_contour(exps(i), ...
                'mark_domain_pts', 0);
        end
    end
    avg_exp.WE_SE.density = avg_A;
    avg_exp.sim_title = exps(1).sim_title;
    % --------------------------------------------------------------------
% %     if plot_all_bouts
% % % %         ALL_bouts = zeros(1e5, 6);
% % % %         ALL_bouts_size = 0;
% % % %         for i = 1:length(exps)
% % % %             % All bouts
% % % %             end_idx = min([length(exps(i).bouts.wake.no_trans.data) ...
% % % %                 length(exps(i).bouts.sleep.no_trans.data)]);
% % % %             bouts_matrix = [exps(i).bouts.wake.no_trans.data(1:end_idx), ...
% % % %                 exps(i).bouts.wake.prev_trans.data(1:end_idx), ...
% % % %                 exps(i).bouts.wake.half_trans.data(1:end_idx), ...
% % % %                 exps(i).bouts.sleep.no_trans.data(1:end_idx), ...
% % % %                 exps(i).bouts.sleep.prev_trans.data(1:end_idx), ...
% % % %                 exps(i).bouts.sleep.half_trans.data(1:end_idx)];
% % % %             nr_bouts = length(bouts_matrix);
% % % %             if ALL_bouts_size + nr_bouts <= length(ALL_bouts)
% % % %                 ALL_bouts(ALL_bouts_size+1:ALL_bouts_size+nr_bouts, :) = bouts_matrix;
% % % %                 ALL_bouts_size = ALL_bouts_size + nr_bouts;
% % % %             end
% % % %         end
% % % %         ALL_bouts(ALL_bouts_size+1:end,:)=[];
% % % %         
% % % %         avg_exp.bouts.wake.no_trans.data = ALL_bouts(:,1);
% % % %         avg_exp.bouts.wake.prev_trans.data = ALL_bouts(:,2);
% % % %         avg_exp.bouts.wake.half_trans.data = ALL_bouts(:,3);
% % % %         avg_exp.bouts.sleep.no_trans.data = ALL_bouts(:,4);
% % % %         avg_exp.bouts.sleep.prev_trans.data = ALL_bouts(:,5);
% % % %         avg_exp.bouts.sleep.half_trans.data = ALL_bouts(:,6);
% %     end
    % --------------------------------------------------------------------
    % Plot average density
    if plot_avg_density
        plot_activity_domains_on_surface_and_contour(avg_exp, ...
            'mark_domain_pts', 0);
    end
    % --------------------------------------------------------------------
    if plot_individual_PL_results
        plot_sleep_wake_bouts_PL_results(exps);
    end
    % --------------------------------------------------------------------
    if display_mean_PL_results
        avg_exp = display_sleep_wake_bouts_PL_results(exps, ...
            display_individual_PL_results, avg_exp);
        avg_exp = display_sleep_wake_bouts_exponential_tail_results(exps, ...
            display_individual_exp_tail_results, avg_exp);
    end
    % --------------------------------------------------------------------
    % Plot bouts (currently only half bouts)
    if plot_all_bouts
        avg_exp = collect_all_bouts_in_avg_exp(exps, avg_exp);
        plot_sleep_wake_bouts_PL_results(avg_exp, 0);
        plot_sleep_wake_bouts_PL_results(avg_exp, 1);
        plot_mean_sleep_wake_bout_results_on_single_figure(exps);
    end
end
% ------------------------------------------------------------------------
end
% folder_name = 'saved_workspaces/fixed_graphs/';
% filenames{5}= [folder_name 'mult_network_sample_105_50_sims_lambda_ratio_20_Mar12_090303.mat'];
% filenames{4}= [folder_name 'mult_network_sample_103_50_sims_lambda_ratio_20_Mar13_175642.mat'];
% filenames{3}= [folder_name 'mult_network_sample_103_50_sims_lambda_ratio_20_Mar13_175642.mat'];
% filenames{2}= [folder_name 'mult_network_sample_102_50_sims_lambda_ratio_20_Mar14_035000.mat'];
% filenames{1}= [folder_name 'mult_network_sample_101_50_sims_lambda_ratio_20_Mar15_115426.mat'];

% folder_name = 'saved_workspaces/fixed_graphs/';
% filenames{17}= [folder_name 'mult_network_sample_5_50_sims_lambda_ratio_20_Mar13_091859.mat'];
% filenames{16}= [folder_name 'mult_network_sample_4_10_sims_lambda_ratio_20_Mar09_212841.mat'];
% filenames{15}= [folder_name 'mult_network_sample_105_10_sims_lambda_ratio_20_Mar03_182005.mat'];
% filenames{14}= [folder_name 'mult_network_sample_104_10_sims_lambda_ratio_20_Mar03_211101.mat'];
% filenames{13}= [folder_name 'mult_network_sample_103_10_sims_lambda_ratio_20_Mar06_225409.mat'];
% filenames{12}= [folder_name 'mult_network_sample_102_10_sims_lambda_ratio_20_Mar04_034216.mat'];
% filenames{11}= [folder_name 'mult_network_sample_101_10_sims_lambda_ratio_20_Mar04_205713.mat'];
% filenames{10}= [folder_name 'mult_network_sample_10_10_sims_lambda_ratio_20_Mar02_204453.mat'];
% filenames{5}= [folder_name 'mult_network_sample_9_10_sims_lambda_ratio_20_Mar02_203222.mat'];
% filenames{9}= [folder_name 'mult_network_sample_8_10_sims_lambda_ratio_20_Feb27_233008.mat'];
% filenames{4}= [folder_name 'mult_network_sample_7_10_sims_lambda_ratio_20_Mar03_214918.mat'];
% filenames{8}= [folder_name 'mult_network_sample_6_10_sims_lambda_ratio_20_Feb28_022247.mat'];
% filenames{3}= [folder_name 'mult_network_sample_5_10_sims_lambda_ratio_20_Mar06_184422.mat'];
% filenames{7}= [folder_name 'mult_network_sample_4_10_sims_lambda_ratio_20_Feb26_173637.mat'];
% filenames{2}= [folder_name 'mult_network_sample_3_10_sims_lambda_ratio_20_Feb26_221036.mat'];
% filenames{6}= [folder_name 'mult_network_sample_2_10_sims_lambda_ratio_20_Mar01_083732.mat'];
% filenames{1}= [folder_name 'mult_network_sample_1_10_sims_lambda_ratio_20_Mar01_192200.mat'];
% 
% folder_name = 'saved_workspaces/random_graphs/';
% filenames{15}= [folder_name 'mult_network_two-state_10_sims_ER_100_9.58_14_SF_100_9.58_14_lambda_ratio_20_Mar02_185908.mat'];
% filenames{14}= [folder_name 'mult_network_two-state_10_sims_ER_100_7.68_12_SF_100_7.68_12_lambda_ratio_20_Mar02_221249.mat'];
% filenames{13}= [folder_name 'mult_network_two-state_10_sims_ER_100_5.78_10_SF_100_5.78_10_lambda_ratio_20_Mar03_022050.mat'];
% filenames{12}= [folder_name 'mult_network_two-state_10_sims_ER_100_3.88_8_SF_100_3.88_8_lambda_ratio_20_Mar03_040813.mat'];
% filenames{11}= [folder_name 'mult_network_two-state_10_sims_ER_100_1.98_6_SF_100_1.98_6_lambda_ratio_20_Mar03_081209.mat'];
% filenames{10}= [folder_name 'mult_network_two-state_10_sims_SF_100_9.58_14_SF_100_9.58_14_lambda_ratio_20_Mar01_010928.mat'];
% filenames{9}= [folder_name 'mult_network_two-state_10_sims_SF_100_7.68_12_SF_100_7.68_12_lambda_ratio_20_Feb27_032701.mat'];
% filenames{8}= [folder_name 'mult_network_two-state_10_sims_SF_100_5.78_10_SF_100_5.78_10_lambda_ratio_20_Feb27_041503.mat'];
% filenames{7}= [folder_name 'mult_network_two-state_10_sims_SF_100_3.88_8_SF_100_3.88_8_lambda_ratio_20_Feb25_051202.mat'];
% filenames{6}= [folder_name 'mult_network_two-state_10_sims_SF_100_1.98_6_SF_100_1.98_6_lambda_ratio_20_Feb27_010148.mat'];
% filenames{5}= [folder_name 'mult_network_two-state_10_sims_ER_100_9.58_14_ER_100_9.58_14_lambda_ratio_20_Feb28_203133.mat'];
% filenames{4}= [folder_name 'mult_network_two-state_10_sims_ER_100_7.68_12_ER_100_7.68_12_lambda_ratio_20_Mar05_131452.mat'];
% filenames{3}= [folder_name 'mult_network_two-state_10_sims_ER_100_5.78_10_ER_100_5.78_10_lambda_ratio_20_Feb27_042832.mat'];
% filenames{2}= [folder_name 'mult_network_two-state_10_sims_ER_100_3.88_8_ER_100_3.88_8_lambda_ratio_20_Feb25_070637.mat'];
% filenames{1}= [folder_name 'mult_network_two-state_10_sims_ER_100_1.98_6_ER_100_1.98_6_lambda_ratio_20_Feb27_094157.mat'];
