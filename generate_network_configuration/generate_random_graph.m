function [mean_degree, adjacency_matrix] = generate_random_graph(varargin)
%% Input parameters
type = varargin{1};
nr_nodes = varargin{2};
par1 = varargin{3};
display_all_results = 0;
mean_degree_tolerance = 0.015;

i=4;
while i<=length(varargin),
    switch varargin{i},
        case 'display_all_results'
            display_all_results = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

display_best_mean_degree_so_far_results = 1;
fprintf('Generating %s random graph with parameters...\n', type);

tSIM = tic;
N = nr_nodes;
adjacency_matrix = zeros(N);

%% Model
if strcmp(type,'ER') 
    %% Erdos-Renyi
    if par1 > 1
        par1 = par1/(nr_nodes-1);
    end
    par1_hat = 2*round(N*(N-1)/2*par1)/(N*(N-1));
    nr_edges_hat = round(N*(N-1)/2*par1);
    mean_degree_hat = 2*round(N*(N-1)/2*par1)/N;
    
    fprintf('\tN = %i, #edges = %i, p = %1.3f, d = %2.3f \n', ...
        N, nr_edges_hat, par1_hat, mean_degree_hat);
    
    if nr_edges_hat ~= N-1
        p_N = par1_hat;
        flag = true;
        current_best_mean_degree_error = 1;
        
        counter = 0;
        while flag
            counter = counter+1;
            temp = rand(N)<p_N;
            temp = triu(temp,1);
            adjacency_matrix = temp + temp';
            mean_degree_error = abs(mean(sum(adjacency_matrix),2)-mean_degree_hat);
            
            % Check the adjacency matrix against connectedness, mean degree
            % tolerance and all-non-zero degrees
            if ~isempty(find(~sum(adjacency_matrix), 1))
                flag = true;
                if display_all_results
                    disp('Zero degree!');
                end
            elseif mean_degree_error ~= 0
                flag=true;
                % Report best mean degree found, since this consumes a good
                % amount of simulation time
                if current_best_mean_degree_error > mean_degree_error
                    if isgraphconnected(adjacency_matrix)
                        if display_all_results
                            disp('Not within better mean degree range!');
                        end
                        if display_best_mean_degree_so_far_results
                            current_best_mean_deg = ...
                                mean(sum(adjacency_matrix), 2);
                            current_best_mean_degree_error = ...
                                mean_degree_error;
                            fprintf('#%i:\t', counter);
                            fprintf('Best mean degree so far is ');
                            fprintf('%2.2f ', current_best_mean_deg);
                            fprintf('with %i edges\n', ...
                                sum(sum(adjacency_matrix))/2);
                        end
                    end
                end
            else
                flag = ~isgraphconnected(adjacency_matrix);
                if display_all_results
                    disp('Checking for connectedness!!!');
                end
            end
        end
    else
        % USE WILSON'S LOOP-ERASING RANDOM WALK APPROACH TO GENERATE A
        % UNIFORM SPANNING TREE OF THE COMPLETE GRAPH 
        fprintf('Using Wilson''s loop-erasing random walk approach\n');
        fprintf('to generate a uniform spanning tree of the complete\n');
        fprintf('graph. This works faster when the graph is a tree.\n');
        
        adjacency_matrix = zeros(N);
        discovered_nodes = zeros(N,1);
        
        nr_discovered_nodes = 1;
        current_node = 1;
        discovered_nodes(nr_discovered_nodes) = 1;
        while nr_discovered_nodes < N
            next_node_candidate_set = setdiff(1:N, current_node);
            next_node_idx = randi(N-1);
            next_node = next_node_candidate_set(next_node_idx);
            if ~ismember(next_node, discovered_nodes(1:nr_discovered_nodes))
                % Update discovered nodes with the next node
                nr_discovered_nodes = nr_discovered_nodes+1;
                discovered_nodes(nr_discovered_nodes) = next_node;
                adjacency_matrix(current_node, next_node) = 1;
                adjacency_matrix(next_node, current_node) = 1;
                if display_all_results
                    fprintf('%3i:(%3i,%3i) is formed\n', ...
                        nr_discovered_nodes, current_node, next_node);
                end
            end
            current_node = next_node;
        end
        
        if ~isgraphconnected(adjacency_matrix)
            error('Generated tree is not connected!!!');
        end
    end
    
elseif strcmp(type,'SF')
    %% 
    par1_hat = ceil(par1*0.5);
    % Preferential attachment currently can wire up to 5 edges for every
    % new vertex
    if ~ismember(par1_hat, [1 2 3 4 5])
        error('Scale-free graph can have 2,4,6,8,10 average connections!!!');
    end    
    
    seed =[ 0 1 0 0 1; 1 0 0 1 0; 0 0 0 1 0; 0 1 1 0 0; 1 0 0 0 0];
    mlinks = par1_hat;
    flag = true;
    while flag
        %         % Some code I found online
        %         adjacency_matrix = SFNG(N,mlinks,seed);

        % My own preferential attachment model
        adjacency_matrix = zeros(N);
        adjacency_matrix(1:5,1:5) = seed;
        current_size = 5;
        while current_size ~= N
            d_k = sum(adjacency_matrix(1:current_size,1:current_size));
            q_k = [0 cumsum(d_k)/sum(d_k)];
            nr_added_edges = 0;
            while nr_added_edges < mlinks
                r = rand;
                i0 = find(r<=q_k,1)-1;
                if adjacency_matrix(current_size+1, i0)==0
                    adjacency_matrix(current_size+1, i0) = 1;
                    adjacency_matrix(i0, current_size+1) = 1;
                    nr_added_edges = nr_added_edges + 1;
                    if display_all_results
                        fprintf('(%3i,%3i) is formed\n', current_size+1, i0);
                    end
                end
            end
            current_size = current_size+1;
        end
        flag=~isgraphconnected(adjacency_matrix);
    end
    
    fprintf('BA algorithm may change the mean degree.\n');
    fprintf('Intially = %1.2f, Now = %1.2f\n', ...
        par1, sum(sum(adjacency_matrix))/N);
    
elseif strcmp(type, 'SF-Chung-Lu')
    %% Added on 10-9-2016
    addpath generate_network_configuration/SF_Chung_Lu/
    addpath ../../power-law-estimator/
    
    mean_degree_goal = par1;  
    nr_edges_hat = mean_degree_goal*N;
    if nr_edges_hat~=round(nr_edges_hat)
        error('FO: Check network parameter!');
    end
    % Determine alpha such that discrete EPL1 distribution on [1, n-1] has
    % the same mean as the input
    orig_mean_deg_goal = mean_degree_goal;
    [alpha, mean_degree_goal] = determine_alpha(N, orig_mean_deg_goal, 0);
    fprintf('Mean degree may have been refined.\n');
    fprintf('Initially = %1.4f, Now = %1.4f\n', ...
        orig_mean_deg_goal, mean_degree_goal);
    
    % Modify network parameters accordinly
    floor_val = floor(mean_degree_goal*N);
    ceil_val = ceil(mean_degree_goal*N);
    if mod(floor_val, 2) == 0
        nr_edges_hat = floor_val;
    elseif mod(ceil_val, 2) == 0
        nr_edges_hat = ceil_val;
    end
    mean_degree_goal = nr_edges_hat/N;
    
    fprintf('\tN = %i, #edges = %i, alpha = %1.2f, d = %2.3f \n', ...
        N, nr_edges_hat, alpha, mean_degree_goal);
    
    current_best_mean_degree_error = inf;
    counter = 0;
    while 1
        counter = counter+1;
        [adjacency_matrix, degrees] = ...
            scale_free_graph_generation(N, alpha, 'zero degree is not ok');
        mean_degree_error = abs(mean(degrees)-mean_degree_goal);
        if current_best_mean_degree_error > mean_degree_error
            if isgraphconnected(adjacency_matrix)
                if display_best_mean_degree_so_far_results
                    current_best_mean_deg = mean(degrees);
                    current_best_mean_degree_error = mean_degree_error;
                    fprintf('#%i:\t', counter);
                    fprintf(['Best mean degree so far is %2.2f ' ...
                        'with %i edges\n'], current_best_mean_deg, ...
                        sum(sum(adjacency_matrix))/2);
                end
                if mean_degree_error < mean_degree_tolerance
                    break;
                end
            end
        end
    end   
        
elseif strcmp(type,'Star')
    adjacency_matrix(1,2:end) = 1;
    
elseif strcmp(type,'Special tree')
    for i=1:floor((N-1)/2)
        adjacency_matrix(i,2*i:2*i+1)=1;
    end
    adjacency_matrix(floor(N/2),N)=1;
    
else
    error('WHAT THE HELL IS THE NETWORK???');
    
end

% adjacency_matrix
mean_degree = sum(sum(adjacency_matrix))/N;
% mean degree for Erdos-Renyi graph is N*p(N)

fprintf('\tDesired mean degree achieved = %6.4f\n', mean_degree);
fprintf('\tTime spent = %3.2f minutes\n\n',toc(tSIM)/60);
end

%% Is the graph connnected?
function connected_or_not = isgraphconnected(G)
% G is a (symmetric) adjaceny matrix
N = length(G);
current_comp = 1;
checked=[];
while ~isempty(setdiff(current_comp,checked))
    temp = setdiff(current_comp,checked);
    current_comp = union(find(G(temp(1),:)==1),current_comp);
    checked = union(checked,temp(1));
end
connected_or_not = length(current_comp)==length(G(1,:));
end
