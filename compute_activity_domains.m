function varargout = compute_activity_domains(varargin)
addpath compute_activity_domains/

%% Input arguments
nW = varargin{1}(1);
nS = varargin{1}(2);
nr_events = varargin{2};
sim_title = varargin{3};
t = varargin{4};
WE = varargin{5};
SE = varargin{6};
A = varargin{7};
plot_pdfs_for_all_intervals = varargin{8};
compute_domains = varargin{9};

%% Model
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

%% Output arguments
varargout{1} = duration_matrix;
varargout{2} = wake_active_domain;
varargout{3} = sleep_active_domain;

end