function varargout = fit_PL_and_EXP_to_SW_bouts(varargin)
addpath ../../power-law-estimator/

%% Input arguments
sleep_bouts_prev_trans = varargin{1};
wake_bouts_prev_trans = varargin{2};
sleep_intervals = varargin{3};
wake_intervals = varargin{4};
sleep_bouts_half_trans = varargin{5};
wake_bouts_half_trans = varargin{6};
display_p_val_stuff = varargin{7};
nr_reps = varargin{8};
fit_PL = varargin{9};
fit_exp = varargin{10};
need_apKS = varargin{11};
sim_title = varargin{12};

%% Model

% Fit bounded power-law to bouts
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

% Fit exponential distribution to bout tails
if fit_exp
    exp_tail_sbpt = fit_exponential_to_tail(sleep_bouts_prev_trans, ...
        'plot_stuff', 0, 'display_p_val_stuff', display_p_val_stuff);
    exp_tail_wbpt = fit_exponential_to_tail(wake_bouts_prev_trans, ...
        'plot_stuff', 0, 'display_p_val_stuff', display_p_val_stuff);
    exp_tail_sbnt = fit_exponential_to_tail(sleep_intervals, ...
        'plot_stuff', 0, 'display_p_val_stuff', display_p_val_stuff);
    exp_tail_wbnt = fit_exponential_to_tail(wake_intervals, ...
        'plot_stuff', 0, 'display_p_val_stuff', display_p_val_stuff);
    exp_tail_sbht = fit_exponential_to_tail(sleep_bouts_half_trans, ...
        'plot_stuff', 0, 'display_p_val_stuff', display_p_val_stuff);
    exp_tail_wbht = fit_exponential_to_tail(wake_bouts_half_trans, ...
        'plot_stuff', 0, 'display_p_val_stuff', display_p_val_stuff);
else
    exp_tail_sbpt = [];
    exp_tail_wbpt = [];
    exp_tail_sbnt = [];
    exp_tail_wbnt = [];
    exp_tail_sbht = [];
    exp_tail_wbht = [];
end

%% Plot bounded power-law fits
temp_SW_exp.bouts.wake.no_trans.data = wake_intervals;
temp_SW_exp.bouts.wake.prev_trans.data = wake_bouts_prev_trans;
temp_SW_exp.bouts.wake.half_trans.data = wake_bouts_half_trans;
temp_SW_exp.bouts.sleep.no_trans.data = sleep_intervals;
temp_SW_exp.bouts.sleep.prev_trans.data = sleep_bouts_prev_trans;
temp_SW_exp.bouts.sleep.half_trans.data = sleep_bouts_half_trans;

temp_SW_exp.bouts.wake.no_trans.PL = PL_wbnt;
temp_SW_exp.bouts.wake.prev_trans.PL = PL_wbpt;
temp_SW_exp.bouts.wake.half_trans.PL = PL_wbht;
temp_SW_exp.bouts.sleep.no_trans.PL = PL_sbnt;
temp_SW_exp.bouts.sleep.prev_trans.PL = PL_sbpt;
temp_SW_exp.bouts.sleep.half_trans.PL = PL_sbht;

temp_SW_exp.bouts.wake.prev_trans.exp_tail = exp_tail_wbpt;
temp_SW_exp.bouts.wake.no_trans.exp_tail = exp_tail_wbnt;
temp_SW_exp.bouts.wake.half_trans.exp_tail = exp_tail_wbht;
temp_SW_exp.bouts.sleep.no_trans.exp_tail = exp_tail_sbnt;
temp_SW_exp.bouts.sleep.prev_trans.exp_tail = exp_tail_sbpt;
temp_SW_exp.bouts.sleep.half_trans.exp_tail = exp_tail_sbht;

temp_SW_exp.sim_title = sim_title;

plot_SW_bouts_PL_results(temp_SW_exp);

%% Output arguments
varargout{1} = PL_sbpt;
varargout{2} = PL_wbpt;
varargout{3} = PL_sbnt;
varargout{4} = PL_wbnt;
varargout{5} = PL_sbht;
varargout{6} = PL_wbht;
varargout{7} = exp_tail_sbpt;
varargout{8} = exp_tail_wbpt;
varargout{9} = exp_tail_sbnt;
varargout{10} = exp_tail_wbnt;
varargout{11} = exp_tail_sbht;
varargout{12} = exp_tail_wbht;

end