function varargout = determine_alpha(varargin)
% DETERMINE_ALPHA finds the exponent alpha for a discrete EPL1 on [1, n-1]
% whose mean is given. This will be used to generate a scale-free network
% with a given mean degree for sleep-wake simulations 

%% Input parameters
n = varargin{1};
goal_mean = varargin{2};
display_stuff = varargin{3};

%% Find alpha
alpha_vec = (1:0.01:5)';
na = length(alpha_vec);
mean_vec = zeros(na, 1);
for i=1:na    
    C = (sum((1:n-1).^(-alpha_vec(i))))^(-1);
    mean_vec(i) = C*sum((1:n-1).^(-alpha_vec(i)+1));
end
% Find the minimum the distance from the goal mean
[~, I] = min(abs(mean_vec-goal_mean));
alpha = alpha_vec(I);
mean = mean_vec(I);

%% Display
if display_stuff
    fprintf('         Goal mean: %3.2f\n', goal_mean);
    fprintf('Power-law exponent: %3.2f\n', alpha_vec(I));
    fprintf('    Power-law mean: %3.2f\n', mean_vec(I));
end

%% Outputs
varargout{1} = alpha;
end