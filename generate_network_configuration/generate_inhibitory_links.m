function [mean_inh_deg, inh_mat] = generate_inhibitory_links(varargin)
% ------------------------------------------------------------------------
nA = varargin{1}(1);            % Size of cluster A
nO = varargin{1}(2);            % Size of (O)pposite cluster
inhPar = varargin{2};           % Inhibition parameter
AA = varargin{3};               % Adjacency matrix of the cluster for which 
                                % an inhibition matrix is constructed
OO = varargin{4};               % Adjacency matrix of the (O)pposite cluster
rule = varargin{5};             % Inhibition rule
% mean_tolerance = varargin{6};   % Tolerance on the mean inhibition degree
% ------------------------------------------------------------------------
plot_inhibitory_degree_distributions = 0;
display_best_mean_degree_so_far_results = 1;
% ------------------------------------------------------------------------
inh_mat = zeros(nA, nO);

if inhPar >= 1 % then it's the average number of inhibition degree
    inhPar = inhPar / nO;
end

current_best_mean_degree_error = inf;
if strcmp(rule, 'ER-like')
    % Error is q is not a probability
    if inhPar<0 || inhPar>1
        error('The inhibition probabilit parameter q must be in [0,1].');
    end
    
    nr_inh_edges_hat = round(nA*nO*inhPar);
    mean_inh_deg_hat = nr_inh_edges_hat/nA;
    q_N = nr_inh_edges_hat/(nA*nO);
    
    
    % Give summary of construction
    fprintf('\nGenerating %ix%i inhibition matrix with parameters...\n', ...
        nA, nO);
    fprintf('\tN = %i, #inh edges = %i, q = %1.4f, di = %2.3f\n', ...
        nA, nr_inh_edges_hat, q_N, mean_inh_deg_hat);
    
    while 1
        % Inhibitory matrix generation
        inh_mat = rand(nA,nO) < inhPar;
        mean_inh_degree = mean(sum(inh_mat,2));
        mean_degree_error = abs(mean_inh_degree - mean_inh_deg_hat);
        % Report best mean degree found, since this consumes a good amount
        % of simulation time
        if current_best_mean_degree_error > mean_degree_error && ...
                display_best_mean_degree_so_far_results
            current_best_mean_degree = mean_inh_degree;
            current_best_mean_degree_error = mean_degree_error;
            fprintf('\tBest mean inhibition degree so far = %6.4f\n', ...
                current_best_mean_degree);
        end
        
        if sum(sum(inh_mat)) == nr_inh_edges_hat
            break;
        end
    end
elseif strcmp(rule, 'Scale-free')
    % 1) Construct degree distributions
    % 2) Form inhibitory links uniformly
    
    % Construct power-law distribution for inhibitory link construction
    alpha = 3-1/(nO+1-inhPar);    % exponent of power-law distribution
    inh_deg_CDF = ((1:nO+1)').^(1-alpha);   
    inh_deg_CDF = cumsum(inh_deg_CDF)./sum(inh_deg_CDF); % CDF
    degrees = (0:nO)';
    
    fprintf('\nGenerating %ix%i inhibition matrix with parameters...\n', nA, nO);
    fprintf([rule '\n']);
    fprintf('dI in (%6.4f,%6.4f), with exponent %1.2f\n', ...
        inhPar - mean_tolerance, inhPar + mean_tolerance, alpha);
    
    while 1
        % Construct degrees from above distribution
        random_inh_degree_probs = rand(nA,1);
        random_inh_degrees = zeros(nA,1);
        for i=1:nA
            random_inh_degrees(i) = ...
                degrees(find(inh_deg_CDF>random_inh_degree_probs(i),1));
        end
        
        % Check mean degree
        mean_degree_error = abs(mean(random_inh_degrees)-inhPar);
        
        % Report best mean degree found, since this consumes a good amount
        % of simulation time
        if current_best_mean_degree_error > mean_degree_error
            current_best_mean_degree = mean(random_inh_degrees);
            current_best_mean_degree_error = mean_degree_error;
            fprintf('\tBest mean inhibition degree so far = %6.4f\n', ...
                current_best_mean_degree);
        end
        
        if mean_degree_error < mean_tolerance
            break;
        end
    end
    
    % Construct inhibition matrix
    inh_mat = zeros(nA,nO);
    for i=1:nA
        which_nodes_are_inhibited = randperm(nO,random_inh_degrees(i));
        inh_mat(i,which_nodes_are_inhibited) = 1;
    end
elseif strcmp(rule, 'Excitatory and inhibitory strengths are same for each node')
    deg = sum(AA,2);
    inh_mat = rand(N)< (deg/(N-1) * ones(1,N));
    while abs(mean(sum(inh_mat,2))-mean(deg))>mean_tolerance
        inh_mat = rand(N)< (deg/(N-1) * ones(1,N));
    end
elseif strcmp(rule, 'Sleep is star')
    inh_mat(:,1)=1;
elseif strcmp(rule, 'Wake is star')
    inh_mat(1,:)=1;
end
% ------------------------------------------------------------------------
if plot_inhibitory_degree_distributions
    hist(sum(inh_mat,2));
    title('Inhibitory degree distributions');
end
% ------------------------------------------------------------------------
% mean inhibitory degree
mean_inh_deg = sum(sum(inh_mat))/nA;

end