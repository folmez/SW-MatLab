function varargout = generate_network_configuration(varargin)
addpath generate_network_configuration/

%% Input arguments
use_sample_experiment = varargin{1};

if use_sample_experiment == 0
    wake_graph = varargin{3};
    nW = varargin{4}(1);
    parW = varargin{4}(2);
    inh_parW = varargin{4}(3);
    sleep_graph = varargin{5};
    nS = varargin{6}(1);
    parS = varargin{6}(2);
    inh_parS = varargin{6}(3);
    WS_inh_rule = varargin{7};
    SW_inh_rule = varargin{8};
    use_undirected_inhibitory_edges_accross_graphs = varargin{9};
    nr_events = varargin{10};
    which_process = varargin{11};
    plot_degree_dist = varargin{12};
elseif use_sample_experiment == 1
    which_sample_experiment = varargin{2};
end

%% Model
if use_sample_experiment == 0
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
        if nW*inh_parW == nS*inh_parS
            [dIW, WS] = generate_inhibitory_links([nW nS], inh_parW, ...
                WW, SS, WS_inh_rule);
            SW = WS';
            dIS = sum(sum(SW))/nS;
        else
            error(['Number of wake-to-sleep inhibitory edges ' ...
                '(%ix%1.2f=%i) must be equal to the number of ' ...
                'sleep-to-wake inhibitory edges (%ix%1.2f=%i)!'], ...
                nW, inh_parW, nW*inh_parW, nS, inh_parS, nS*inh_parS);
        end
    else
        [dIW, WS] = generate_inhibitory_links([nW nS], inh_parW, ...
            WW, SS, WS_inh_rule);
        [dIS, SW] = generate_inhibitory_links([nS nW], inh_parS, ...
            SS, WW, SW_inh_rule);
    end
    
    % Set simulation title
    sim_title = ['W=' wake_graph '(' num2str(nW) ',' ...
        num2str(dW, '%2.2f') ',' num2str(dIW, '%2.2f') ...
        '), S=' sleep_graph '(' num2str(nS) ',' num2str(dS, '%2.2f') ...
        ',' num2str(dIS, '%2.2f') ')' ', # of events=' ...
        num2str(nr_events,'%1.0e, ') which_process];
    
elseif use_sample_experiment == 1
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

%% Plot degree distributions
if plot_degree_dist
    figure,

    subplot(2,2,1),
    histogram(sum(WW));
    title('Wake degree dist');
    legend(['de_W=' num2str(dW)]);
    
    subplot(2,2,2),
    histogram(sum(SS));
    title('Sleep degree dist');
    legend(['de_S=' num2str(dS)]);
    
    subplot(2,2,3),
    histogram(sum(WS,2));
    title('W->S inh. deg. dist.');
    legend(['di_W_-_>_S=' num2str(dIW)]);
    
    subplot(2,2,4),
    histogram(sum(SW,2));
    title('S->W inh. deg. dist.');
    legend(['di_S_-_>_W=' num2str(dIS)]);
    
    drawnow;
end

%% Output arguments
varargout{1} = nr_events;
varargout{2} = sim_title;
varargout{3} = WW;
varargout{4} = dW;
varargout{5} = SS;
varargout{6} = dS;
varargout{7} = wake_graph;
varargout{8} = sleep_graph;
varargout{9} = nW;
varargout{10} = nS;
varargout{11} = SW;
varargout{12} = dIS;
varargout{13} = WS;
varargout{14} = dIW;

end