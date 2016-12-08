function [wad, sad] = compute_wake_sleep_domains_NEW(varargin)
%   [wad, sad] = COMPUTE_WAKE_SLEEP_DOMAINS_NEW([nW nS], sim_title, A);
%   COMPUTE_WAKE_SLEEP_DOMAINS_NEW computes wake and sleep domains for a 
%   given simulated (WE(t),SE(t)). An error is displayed if there is only 
%   one state.
%% Input parameters
nW = varargin{1}(1);
nS = varargin{1}(2);
sim_title = varargin{2};
A = varargin{3};

i = 8;
while i<=length(varargin),
    switch varargin{i},
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

% STAGE 1
% Rough detection of domains
nr_test_levels = 50;
levels = unique(A(:));
nr_levels = length(levels);
test_level_idx = 1:round(nr_levels/nr_test_levels):nr_levels;
test_level_connectedness = (-1)*ones(nr_levels, 1);

plot_levels = 0;
for i = test_level_idx
    [x, y, ~] = find(A>=levels(i));
    
    if plot_levels
        figure, contour(A,50); hold,
        plot(y, x, 'kx');
        title(['Level ' num2str(i)], 'FontSize', 20);
        h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
        h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
    end
    
    test_level_connectedness(i) = is_domain_connected(nW, nS, x, y, i);
    
    fprintf('Level %4i/%4i: ', i, nr_levels);
    if test_level_connectedness(i)
        fprintf('Connected\n');
    else
        fprintf('NOT connected\n');
    end
    
    if plot_levels
        pause;
        close;
    end
end
fprintf('\n');

% STAGE 2
% Fine detection domains
not_conn_test_level_idx = find(test_level_connectedness(test_level_idx)==0);
nr_temp = length(not_conn_test_level_idx);
temp = zeros(nr_temp, 1);
for i=1:nr_temp
    j = i;
    while not_conn_test_level_idx(i)+(j-i) == not_conn_test_level_idx(j)
        temp(i) = temp(i)+1;
        j = j+1;
        if j>nr_temp
            break
        end
    end
end
[~, max_idx] = max(temp);
temp_idx = not_conn_test_level_idx(max_idx);

test_level_idx = test_level_idx(temp_idx-1):test_level_idx(temp_idx);
test_level_connectedness = (-1)*ones(nr_levels, 1);
for i = test_level_idx
    [x, y, ~] = find(A>=levels(i));
    
    if plot_levels
        figure, contour(A,50); hold,
        plot(y, x, 'kx');
        title(['Level ' num2str(i)], 'FontSize', 20);
        h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
        h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
    end
    
    test_level_connectedness(i) = is_domain_connected(nW, nS, x, y, i);
    
    fprintf('Level %4i/%4i: ', i, nr_levels);
    if test_level_connectedness(i)
        fprintf('Connected\n');
    else
        fprintf('NOT connected\n');
    end
    
    if plot_levels
        pause;
        close;
    end
end

temp = test_level_connectedness(test_level_idx);
idx = find(~temp, 1, 'first');
best_two_domain_level = test_level_idx(idx);

[x, y, ~] = find(A>=levels(best_two_domain_level));
% figure, contour(A,50); hold,
% plot(y, x, 'kx');
% title(['Level ' num2str(best_two_domain_level)], 'FontSize', 20);
% h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
% h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);

% A wake domain index and a sleep domain index
[~, wd_idx] = max(x);
wp = [x(wd_idx(1)) y(wd_idx(1))];
[~, sd_idx] = max(y);
sp = [x(sd_idx(1)) y(sd_idx(1))];

wad = find_largest_connected_superset(nW, nS, [x y], wp);
sad = find_largest_connected_superset(nW, nS, [x y], sp);

% Plot wake- and sleep- active domain
figure, contour(A,50); hold,
plot(wad(:,2), wad(:,1), 'kx');
plot(sad(:,2), sad(:,1), 'ko');
title(sim_title, 'FontSize', 15);
h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
h_legend = legend('PDF of SE x WE', ...
    'Wake-active domain', 'Sleep-active domain', 'Location','Best');
set(h_legend,'FontSize',15);

% We need to return wad and sad after subtracting a 1 because A was
% generated after a 1 was added to WE(t) and SE(t)
wad = wad - 1;
sad = sad - 1;
end

%% IS DOMAIN CONNECTED?
function flag = is_domain_connected(nW, nS, x, y, i)
plot_stuff = 0;

B = [x y];
[nr_B, ~ ] = size(B); 
lcss = zeros((nW+1)*(nS+1),2); % Largest Connected Superset

% Pick the first element of subset1
lcss(1,:) = [x(1) y(1)];
nr_lcss = 1;
rap = lcss(1,:); % Recently Added Points
while ~isempty(rap)
    % Find the indices of neighbor points of the recently added points
    % Candidate Subset1 Points
    [nr_rap, ~] = size(rap);
    
    %     cs1p = [rap+ones(nr_rap,1)*[-1 0]; rap+ones(nr_rap,1)*[1 0]; ...
    %         rap+ones(nr_rap,1)*[0 1]; rap+ones(nr_rap,1)*[0 -1]];
    %
    
    %     cs1p = [rap+ones(nr_rap,1)*[-1 0]; rap+ones(nr_rap,1)*[1 0]; ...
    %         rap+ones(nr_rap,1)*[0 1]; rap+ones(nr_rap,1)*[0 -1]; ...
    %         rap+ones(nr_rap,1)*[1 1]; rap+ones(nr_rap,1)*[-1 -1]; ...
    %         rap+ones(nr_rap,1)*[-1 1]; rap+ones(nr_rap,1)*[1 -1]];
    
    cs1p = [rap+ones(nr_rap,1)*[-1 0]; rap+ones(nr_rap,1)*[1 0]; ...
        rap+ones(nr_rap,1)*[0 1]; rap+ones(nr_rap,1)*[0 -1]; ...
        rap+ones(nr_rap,1)*[1 1]; rap+ones(nr_rap,1)*[-1 -1]; ...
        rap+ones(nr_rap,1)*[-1 1]; rap+ones(nr_rap,1)*[1 -1]; ...
        rap+ones(nr_rap,1)*[-2 0]; rap+ones(nr_rap,1)*[2 0]; ...
        rap+ones(nr_rap,1)*[0 2]; rap+ones(nr_rap,1)*[0 -2]; ...
        rap+ones(nr_rap,1)*[2 2]; rap+ones(nr_rap,1)*[-2 -2]; ...
        rap+ones(nr_rap,1)*[-2 2]; rap+ones(nr_rap,1)*[2 -2]; ...
        rap+ones(nr_rap,1)*[-2 1]; rap+ones(nr_rap,1)*[2 1]; ...
        rap+ones(nr_rap,1)*[1 2]; rap+ones(nr_rap,1)*[1 -2]; ...
        rap+ones(nr_rap,1)*[2 -1]; rap+ones(nr_rap,1)*[-2 -1]; ...
        rap+ones(nr_rap,1)*[-1 2]; rap+ones(nr_rap,1)*[-1 -1]];
    
    % Remove repetitions
    cs1p = unique(cs1p, 'rows');
    % Remove out-of-grid points
    cs1p(cs1p(:,1)<1,:) = []; cs1p(cs1p(:,1)>nW+1,:) = [];
    cs1p(cs1p(:,2)<1,:) = []; cs1p(cs1p(:,2)>nS+1,:) = [];

    % Subset1 Points to be added are those candidate subset1 points which
    % are in B but not already in subset1
    s1p_tba = setdiff( intersect(cs1p, B, 'rows'), ...
        lcss(1:nr_lcss,:), 'rows');
    [nr_s1p_tba, ~] = size(s1p_tba);
    % Put them in subset1
    lcss(nr_lcss+1:nr_lcss+nr_s1p_tba, :) = s1p_tba;
    % Update Recently Added Points and the Number of Subset1 Points
    rap = s1p_tba;
    nr_lcss = nr_lcss+nr_s1p_tba;
end

flag = nr_lcss==nr_B;

if plot_stuff 
    figure,
    plot(y, x, 'kx');
    axis([0 nS+2 0 nW+2])
    connectedness_title = ['Level ' num2str(i) ' is '];
    if flag
        title([connectedness_title 'connected'], 'FontSize', 20);
    else
        title([connectedness_title 'NOT connected'], 'FontSize', 20);
    end
    h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
    h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
    pause;
    close
end
end

%% Find largest connected superset
function lcss = find_largest_connected_superset(nW, nS, B, single_element)
% Initiate the first subset
lcss = zeros((nW+1)*(nS+1),2); % Largest Connected Superset

% Pick the first element of subset1
lcss(1,:) = single_element;
nr_lcss = 1;
rap = lcss(1,:); % Recently Added Points
while ~isempty(rap)
    % Find the indices of neighbor points of the recently added points
    % Candidate Subset1 Points
    [nr_rap, ~] = size(rap);
    
    cs1p = [rap+ones(nr_rap,1)*[-1 0]; rap+ones(nr_rap,1)*[1 0]; ...
        rap+ones(nr_rap,1)*[0 1]; rap+ones(nr_rap,1)*[0 -1]];

    % Remove repetitions
    cs1p = unique(cs1p, 'rows');
    % Remove out-of-grid points
    cs1p(cs1p(:,1)<1,:) = []; cs1p(cs1p(:,1)>nW+1,:) = [];
    cs1p(cs1p(:,2)<1,:) = []; cs1p(cs1p(:,2)>nS+1,:) = [];
    
    % Subset1 Points to be added are those candidate subset1 points which
    % are in B but not already in subset1
    s1p_tba = setdiff( intersect(cs1p, B, 'rows'), ...
        lcss(1:nr_lcss,:), 'rows');
    [nr_s1p_tba, ~] = size(s1p_tba);
    % Put them in subset1
    lcss(nr_lcss+1:nr_lcss+nr_s1p_tba, :) = s1p_tba;
    % Update Recently Added Points and the Number of Subset1 Points
    rap = s1p_tba;
    nr_lcss = nr_lcss+nr_s1p_tba;
end

% Remove extra zero rows
lcss(lcss(:,1)==0,:) = [];
end