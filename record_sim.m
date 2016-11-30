function varargout = record_sim(varargin)
%% Input arguments
sim_title = varargin{1};
nr_events = varargin{2};
wake_graph = varargin{3};
nW = varargin{4};
WW = varargin{5};
dW = varargin{6};
WS = varargin{7};
dIW = varargin{8};
sleep_graph = varargin{9};
nS = varargin{10};
SS = varargin{11};
dS = varargin{12};
SW = varargin{13};
dIS = varargin{14};
A = varargin{15};
duration_matrix = varargin{16};
wake_active_domain = varargin{17};
sleep_active_domain = varargin{18};
wake_intervals = varargin{19};
sleep_intervals = varargin{20};
transition_intervals = varargin{21};
wake_bouts_prev_trans = varargin{22};
wake_bouts_half_trans = varargin{23};
sleep_bouts_prev_trans = varargin{24};
sleep_bouts_half_trans = varargin{25};
PL_wbnt = varargin{26};
PL_wbpt = varargin{27};
PL_wbht = varargin{28};
PL_sbnt = varargin{29};
PL_sbpt = varargin{30};
PL_sbht = varargin{31};
exp_tail_wbpt = varargin{32};
exp_tail_wbnt = varargin{33};
exp_tail_wbht = varargin{34};
exp_tail_sbnt = varargin{35};
exp_tail_sbpt = varargin{36};
exp_tail_sbht = varargin{37};
F = varargin{38};
D = varargin{39};
F_WE = varargin{40};
D_WE = varargin{41};
F_SE = varargin{42};
D_SE = varargin{43};
use_gillespie_algorithm = varargin{44};

%% Model
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

SW_exp.duration_matrix = duration_matrix;
SW_exp.WE_SE.wake_active_domain = wake_active_domain;
SW_exp.WE_SE.sleep_active_domain = sleep_active_domain;

SW_exp.intervals.wake = wake_intervals;
SW_exp.intervals.sleep = sleep_intervals;
SW_exp.intervals.transition = transition_intervals;

SW_exp.bouts.wake.no_trans.data = wake_intervals;
SW_exp.bouts.wake.prev_trans.data = wake_bouts_prev_trans;
SW_exp.bouts.wake.half_trans.data = wake_bouts_half_trans;
SW_exp.bouts.sleep.no_trans.data = sleep_intervals;
SW_exp.bouts.sleep.prev_trans.data = sleep_bouts_prev_trans;
SW_exp.bouts.sleep.half_trans.data = sleep_bouts_half_trans;

SW_exp.bouts.wake.no_trans.PL = PL_wbnt;
SW_exp.bouts.wake.prev_trans.PL = PL_wbpt;
SW_exp.bouts.wake.half_trans.PL = PL_wbht;
SW_exp.bouts.sleep.no_trans.PL = PL_sbnt;
SW_exp.bouts.sleep.prev_trans.PL = PL_sbpt;
SW_exp.bouts.sleep.half_trans.PL = PL_sbht;

SW_exp.bouts.wake.prev_trans.exp_tail = exp_tail_wbpt;
SW_exp.bouts.wake.no_trans.exp_tail = exp_tail_wbnt;
SW_exp.bouts.wake.half_trans.exp_tail = exp_tail_wbht;
SW_exp.bouts.sleep.no_trans.exp_tail = exp_tail_sbnt;
SW_exp.bouts.sleep.prev_trans.exp_tail = exp_tail_sbpt;
SW_exp.bouts.sleep.half_trans.exp_tail = exp_tail_sbht;

SW_exp.DDE.F = F;
SW_exp.DDE.D = D;
SW_exp.DDE.F_WE = F_WE;
SW_exp.DDE.D_WE = D_WE;
SW_exp.DDE.F_SE = F_SE;
SW_exp.DDE.D_SE = D_SE;

SW_exp.use_gillespie_algorithm = use_gillespie_algorithm;

%% Output arguments
varargout{1} = SW_exp;

end