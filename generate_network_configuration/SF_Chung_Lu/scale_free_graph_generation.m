function varargout = scale_free_graph_generation(varargin)
% Note: Original Chung-Lu model (chung_lu.m) generates 0-degree nodes. This
% function generates a random network with no 0-degree nodes.

addpath ../../power-law-estimator/

%% Model parameters
n = varargin{1};
alpha = varargin{2};
zero_degree_status = varargin{3};   % 'zero degree is ok'
                                    % 'zero degree is not ok'
display_stuff = 0;
plot_stuff = 0;

i = 4;
while i<=length(varargin),
    switch varargin{i},
        case 'plot_stuff',      plot_stuff = varargin{i+1};
        case 'display_stuff',   display_stuff = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

use_discrete_epl1 = 1;

%% Model
new_n = n;
i = 1;
nr_zero_deg_nodes_vec = zeros(n, 1);
if display_stuff
    fprintf('\tSim#\tSize\t#0-deg nodes\n');
end
while 1
    % Generate power-law distributed integers
    while 1
        if use_discrete_epl1
            k = generate_synthetic_data_from_CDF(...
                'EPL1-Discrete', alpha, [1 n-1], new_n, 0);
        else
            k = 0.5*(1-rand(new_n,1)).^(-1/(alpha-1));
        end
        if max(k)<n
            break;
        end
    end
    k = round(k);
    
    % Run Chung-Lu model to generate a random graph based on k
    [A, nr_zero_deg_nodes] = chung_lu(k);
    if strcmp(zero_degree_status, 'zero degree is not ok')
        % Display progress
        if display_stuff
            fprintf('\t%i\t%i\t%i\n', i, new_n, nr_zero_deg_nodes);
        end
        
        % Break if number of number of non 0-degrees nodes is equal to n
        if new_n-nr_zero_deg_nodes == n
            break;
        else
            % Update nr_zero degree nodes vector
            nr_zero_deg_nodes_vec(i) = nr_zero_deg_nodes;
            nr_nodes_to_add = round(mean(nr_zero_deg_nodes_vec(1:i)));
            new_n = n + nr_nodes_to_add;
            i = i+1;
        end
    elseif strcmp(zero_degree_status, 'zero degree is ok')
        break;
    end
end

% Permute adjacency matrix because hubs are at the bottom
node_perm = randperm(n);
A_perm = zeros(n);
for i = 1:n
    A_perm(i,:) = A(node_perm(i), node_perm);
end
A = A_perm;

% Find degrees
degrees = sum(A,2);
if strcmp(zero_degree_status, 'zero degree is not ok')
    % Remove 0-degree rows and columns from the adjacaceny matrix
    A(degrees==0, :) = [];
    A(:, degrees==0) = [];
    degrees(degrees==0) = [];
end

%% Plot stuff
if plot_stuff    
    idx = degrees~=0;
    data_title = ['Avg deg. of ' num2str(n) ' nodes = ' ...
        num2str(mean(degrees)) ...
        ' (' num2str(n-sum(idx)) ' 0-deg nodes)'];
    plot_approximate_pdf_of_data(degrees(idx), ...
        'data_title', data_title, ...
        'power-law fit', [alpha 1 max(degrees) 0]);
end

%% Outputs
varargout{1} = A;
varargout{2} = degrees;

end