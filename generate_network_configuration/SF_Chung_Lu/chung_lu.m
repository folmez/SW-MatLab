%% See https://e-reports-ext.llnl.gov/pdf/802505.pdf for more details
function varargout = chung_lu(varargin)
%% Model parameters
k = varargin{1};    % nx1 vector holding expected degrees

m = sum(k)/2;
n = length(k);      % Number of nodes
plot_stuff = 0;
display_stuff = 0;

%% Check inputs
if max(k)>n
    error('FO: Inputted degrees cannot be larger than n-1!')
end

%% Model
% An edge (i,j)=(j,i) will be formed with probability (k_i*k_j)/(2m). Self
% edges and multiple edges are not allowed.

% Generate a metrix whose entries are above probabilities
k = reshape(k,n,1);
k = sort(k);
B = zeros(n);
for i=1:n
    B(i,i+1:n) = k(i)*k(i+1:n)/(2*m);
end

% Generate a matrix whose entries are from U(0,1).
A = rand(n);
% Consider only the upper triangular section
A = triu(A, 1);

% Turn A into an adjaceny matrix by comparing it against B
A = A<B;

% Sum with its transpose to guarantee symmetricity
A = A + A';

% Calculate simulation degrees
degrees = sort(sum(A))';

%% Display top three hubs
i = 10;
top_inputs = k(end-i+1:end);
top_expecteds = round(k(end-i+1:end) - k(end-i+1:end).^2/(2*m));
top_sims = degrees(end-i+1:end);
nr_zero_deg_nodes = sum(degrees==0);
if display_stuff
    fprintf('\t\tTOP THREE HUBS\n');
    fprintf('\t#\tInput\tExp''ted\tSimulation\n');
    fprintf('\t%i\t%3i\t%3i\t%3i\n', [(i:-1:1)', top_inputs, ...
        top_expecteds, top_sims]');
    fprintf('\n');
    fprintf('       Relative number of\n\t  0-degree nodes:');
    fprintf(' %i/%i = %2.2f', nr_zero_deg_nodes, n, nr_zero_deg_nodes/n);
    fprintf('\n');    
    fprintf('     Average input degree: %3.2f\n', mean(k));
    fprintf('Average simulation degree: %3.2f\n', mean(degrees));
end

%% Plot results
if plot_stuff
    figure,
    %     subplot(2,1,1);
    plot(1:n, degrees - k, 'bo');
    title('Difference between Expected Degree and Degree ');
    %     subplot(2,1,2);
    %     plot(1:n, abs(degrees./k - 1), 'bo');
    %     title('Relative difference between Expected Degree and Degree ');
    xlabel('Node (degree increasing from left to right)');
end

%% Outputs
varargout{1} = A;
varargout{2} = nr_zero_deg_nodes;
end