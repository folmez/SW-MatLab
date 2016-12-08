function varargout = compute_wake_sleep_domains_3(varargin)
% This method, as opposed to the previous activity domain computation
% methods, does not rely on contor lines being connected or anything like
% that. It searches the heat map from bottom to top and splits the region
% into two clusters using MatLab's built-in k-means function. If clusters
% are distant enough, method stops and identifies activity domains.

%% Input parameters
nW = varargin{1}(1);
nS = varargin{1}(2);
sim_title = varargin{2};
A = varargin{3};

plot_all_results = 0;
plot_activity_domains = 1;

%% Step 1: Find activity domains

% Find all levels of the heat map
levels = unique(A(:));
% Start a figure to see the progress of the activity domains
if plot_all_results, figure, end
for i = 1:length(levels)
    % Find current level
    current_level = levels(i);
    % Find coordinates that are not below the current level
    [x, y, ~] = find( A >= current_level );
    % Split into two clusters if there is more than one point
    if length(x)>1
        % Split using k-means
        [idx, C] = ...
            kmeans([y x], 2, 'Distance', 'cityblock','Replicates', 5);
        
        % Plot clusters
        if plot_all_results
            plot(y(idx==1), x(idx==1), 'r.', 'MarkerSize', 12);
            hold on;
            plot(y(idx==2), x(idx==2), 'b.', 'MarkerSize', 12);
            plot(C(:,1), C(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
            legend('Cluster 1', 'Cluster 2', 'Centroids', ...
                'Location', 'NE');
            xlim([0 nS]);
            ylim([0 nW]);
        end
        
        % Enlarge the clusters by one neighbor to see if they intersect. If
        % they do then we must go up by one more level. If they don't then
        % we have found the activity domains
        x1 = x(idx==1);
        y1 = y(idx==1);
        x2 = x(idx==2);
        y2 = y(idx==2);
        x1_en = [ x1 ; x1-1 ; x1+1 ; x1   ; x1   ];
        y1_en = [ y1 ; y1   ; y1   ; y1-1 ; y1+1 ];
        x2_en = [ x2 ; x2-1 ; x2+1 ; x2   ; x2   ];
        y2_en = [ y2 ; y2   ; y2   ; y2-1 ; y2+1 ];
        
        % Plot title
        if plot_all_results
            if ~isempty(intersect([x1_en y1_en], [x2_en y2_en], 'rows'))
                title(['Level ' num2str(i) ' not split'], 'FontSize', 20);
            else
                title(['Level ' num2str(i) ' SPLIT !!!'], 'FontSize', 20);
            end
            hold off;
            pause;
        end
        
        % Break out of the for loop if the activity domains are found
        if isempty(intersect([x1_en y1_en], [x2_en y2_en], 'rows'))
            break;
        end
    end
end

% Plot axis labels
if plot_all_results
    h_xlabel = xlabel('SE');
    set(h_xlabel,'FontSize',15);
    h_ylabel = ylabel('WE');
    set(h_ylabel,'FontSize',15);
end

% Identify activity domains
if max(x1) > max(x2)
    wad = [x1 y1];
    sad = [x2 y2];
else
    wad = [x2 y2];
    sad = [x1 y1];
end

%% Step 2: Draw a convex hull for the activity domains and consider all points inside the hull

% Define grid
x_mat = (0:nS)'*ones(1,nW+1);
y_mat = ones(nS+1,1)*(0:nW);

% Find convex hulls containing the activity domains
wk = convhull(wad(:,1), wad(:,2));
sk = convhull(sad(:,1), sad(:,2));

% Find insides
in_wad = inpolygon(x_mat, y_mat, wad(wk,2), wad(wk,1));
in_sad = inpolygon(x_mat, y_mat, sad(sk,2), sad(sk,1));

% Plot activity domains
if plot_activity_domains
    figure,   
    plot(wad(:,2), wad(:,1), 'kx');
    hold on;
    plot(wad(wk,2), wad(wk,1), 'r-');
    plot(x_mat(in_wad), y_mat(in_wad), 'k.', 'MarkerSize', 1);    
    plot(sad(:,2), sad(:,1), 'ko');
    plot(sad(sk,2), sad(sk,1), 'r-');
    plot(x_mat(in_sad), y_mat(in_sad), 'k.', 'MarkerSize', 1);        
    title(sim_title, 'FontSize', 15);
    h_xlabel = xlabel('SE');
    set(h_xlabel,'FontSize',15);
    h_ylabel = ylabel('WE');
    set(h_ylabel,'FontSize',15);
    h_legend = legend('Wake-active domain', 'Sleep-active domain', ...
        'Location','Best');
    set(h_legend,'FontSize',15);
    xlim([0 nS]);
    ylim([0 nW]);
end

% Append insides to the activity domains
wad = [y_mat(in_wad) x_mat(in_wad)];
sad = [x_mat(in_wad) y_mat(in_wad)];

if plot_all_results
    figure,
    plot(wad(:,2), wad(:,1), 'kx');
    hold on;
    plot(sad(:,2), sad(:,1), 'ko');
    title(sim_title, 'FontSize', 15);
    h_xlabel = xlabel('SE');
    set(h_xlabel,'FontSize',15);
    h_ylabel = ylabel('WE');
    set(h_ylabel,'FontSize',15);
    h_legend = legend('Wake-active domain', 'Sleep-active domain', ...
        'Location','Best');
    set(h_legend,'FontSize',15);
    xlim([0 nS]);
    ylim([0 nW]);
end

% We need to return wad and sad after subtracting 1 because A was generated
% after 1 was added to WE(t) and SE(t)
wad = wad - 1;
sad = sad - 1;

%% Output arguments
varargout{1} = wad;
varargout{2} = sad;

end

% function [wad, sad] = compute_wake_sleep_domains_NEW(varargin)
% %   [wad, sad] = COMPUTE_WAKE_SLEEP_DOMAINS_NEW([nW nS], sim_title, A);
% %   COMPUTE_WAKE_SLEEP_DOMAINS_NEW computes wake and sleep domains for a 
% %   given simulated (WE(t),SE(t)). An error is displayed if there is only 
% %   one state.
% %% Input parameters
% nW = varargin{1}(1);
% nS = varargin{1}(2);
% sim_title = varargin{2};
% A = varargin{3};

% % STAGE 1
% % Rough detection of domains
% nr_test_levels = 50;
% levels = unique(A(:));
% nr_levels = length(levels);
% test_level_idx = 1:round(nr_levels/nr_test_levels):nr_levels;
% test_level_connectedness = (-1)*ones(nr_levels, 1);
% 
% plot_levels = 0;
% for i = test_level_idx
%     [x, y, ~] = find(A>=levels(i));
%     
%     if plot_levels
%         figure, contour(A,50); hold,
%         plot(y, x, 'kx');
%         title(['Level ' num2str(i)], 'FontSize', 20);
%         h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
%         h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
%     end
%     
%     test_level_connectedness(i) = is_domain_connected(nW, nS, x, y, i);
%     
%     fprintf('Level %4i/%4i: ', i, nr_levels);
%     if test_level_connectedness(i)
%         fprintf('Connected\n');
%     else
%         fprintf('NOT connected\n');
%     end
%     
%     if plot_levels
%         pause;
%         close;
%     end
% end
% fprintf('\n');
% 
% % STAGE 2
% % Fine detection domains
% not_conn_test_level_idx = find(test_level_connectedness(test_level_idx)==0);
% nr_temp = length(not_conn_test_level_idx);
% temp = zeros(nr_temp, 1);
% for i=1:nr_temp
%     j = i;
%     while not_conn_test_level_idx(i)+(j-i) == not_conn_test_level_idx(j)
%         temp(i) = temp(i)+1;
%         j = j+1;
%         if j>nr_temp
%             break
%         end
%     end
% end
% [~, max_idx] = max(temp);
% temp_idx = not_conn_test_level_idx(max_idx);
% 
% test_level_idx = test_level_idx(temp_idx-1):test_level_idx(temp_idx);
% test_level_connectedness = (-1)*ones(nr_levels, 1);
% for i = test_level_idx
%     [x, y, ~] = find(A>=levels(i));
%     
%     if plot_levels
%         figure, contour(A,50); hold,
%         plot(y, x, 'kx');
%         title(['Level ' num2str(i)], 'FontSize', 20);
%         h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
%         h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
%     end
%     
%     test_level_connectedness(i) = is_domain_connected(nW, nS, x, y, i);
%     
%     fprintf('Level %4i/%4i: ', i, nr_levels);
%     if test_level_connectedness(i)
%         fprintf('Connected\n');
%     else
%         fprintf('NOT connected\n');
%     end
%     
%     if plot_levels
%         pause;
%         close;
%     end
% end
% 
% temp = test_level_connectedness(test_level_idx);
% idx = find(~temp, 1, 'first');
% best_two_domain_level = test_level_idx(idx);
% 
% [x, y, ~] = find(A>=levels(best_two_domain_level));
% % figure, contour(A,50); hold,
% % plot(y, x, 'kx');
% % title(['Level ' num2str(best_two_domain_level)], 'FontSize', 20);
% % h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
% % h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
% 
% % A wake domain index and a sleep domain index
% [~, wd_idx] = max(x);
% wp = [x(wd_idx(1)) y(wd_idx(1))];
% [~, sd_idx] = max(y);
% sp = [x(sd_idx(1)) y(sd_idx(1))];
% 
% wad = find_largest_connected_superset(nW, nS, [x y], wp);
% sad = find_largest_connected_superset(nW, nS, [x y], sp);
% 
% % Plot wake- and sleep- active domain
% figure, contour(A,50); hold,
% plot(wad(:,2), wad(:,1), 'kx');
% plot(sad(:,2), sad(:,1), 'ko');
% title(sim_title, 'FontSize', 15);
% h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
% h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
% h_legend = legend('PDF of SE x WE', ...
%     'Wake-active domain', 'Sleep-active domain', 'Location','Best');
% set(h_legend,'FontSize',15);
% 
% % We need to return wad and sad after subtracting a 1 because A was
% % generated after a 1 was added to WE(t) and SE(t)
% wad = wad - 1;
% sad = sad - 1;
% end
% 
% %% IS DOMAIN CONNECTED?
% function flag = is_domain_connected(nW, nS, x, y, i)
% plot_stuff = 0;
% 
% B = [x y];
% [nr_B, ~ ] = size(B); 
% lcss = zeros((nW+1)*(nS+1),2); % Largest Connected Superset
% 
% % Pick the first element of subset1
% lcss(1,:) = [x(1) y(1)];
% nr_lcss = 1;
% rap = lcss(1,:); % Recently Added Points
% while ~isempty(rap)
%     % Find the indices of neighbor points of the recently added points
%     % Candidate Subset1 Points
%     [nr_rap, ~] = size(rap);
%     
%     %     cs1p = [rap+ones(nr_rap,1)*[-1 0]; rap+ones(nr_rap,1)*[1 0]; ...
%     %         rap+ones(nr_rap,1)*[0 1]; rap+ones(nr_rap,1)*[0 -1]];
%     %
%     
%     %     cs1p = [rap+ones(nr_rap,1)*[-1 0]; rap+ones(nr_rap,1)*[1 0]; ...
%     %         rap+ones(nr_rap,1)*[0 1]; rap+ones(nr_rap,1)*[0 -1]; ...
%     %         rap+ones(nr_rap,1)*[1 1]; rap+ones(nr_rap,1)*[-1 -1]; ...
%     %         rap+ones(nr_rap,1)*[-1 1]; rap+ones(nr_rap,1)*[1 -1]];
%     
%     cs1p = [rap+ones(nr_rap,1)*[-1 0]; rap+ones(nr_rap,1)*[1 0]; ...
%         rap+ones(nr_rap,1)*[0 1]; rap+ones(nr_rap,1)*[0 -1]; ...
%         rap+ones(nr_rap,1)*[1 1]; rap+ones(nr_rap,1)*[-1 -1]; ...
%         rap+ones(nr_rap,1)*[-1 1]; rap+ones(nr_rap,1)*[1 -1]; ...
%         rap+ones(nr_rap,1)*[-2 0]; rap+ones(nr_rap,1)*[2 0]; ...
%         rap+ones(nr_rap,1)*[0 2]; rap+ones(nr_rap,1)*[0 -2]; ...
%         rap+ones(nr_rap,1)*[2 2]; rap+ones(nr_rap,1)*[-2 -2]; ...
%         rap+ones(nr_rap,1)*[-2 2]; rap+ones(nr_rap,1)*[2 -2]; ...
%         rap+ones(nr_rap,1)*[-2 1]; rap+ones(nr_rap,1)*[2 1]; ...
%         rap+ones(nr_rap,1)*[1 2]; rap+ones(nr_rap,1)*[1 -2]; ...
%         rap+ones(nr_rap,1)*[2 -1]; rap+ones(nr_rap,1)*[-2 -1]; ...
%         rap+ones(nr_rap,1)*[-1 2]; rap+ones(nr_rap,1)*[-1 -1]];
%     
%     % Remove repetitions
%     cs1p = unique(cs1p, 'rows');
%     % Remove out-of-grid points
%     cs1p(cs1p(:,1)<1,:) = []; cs1p(cs1p(:,1)>nW+1,:) = [];
%     cs1p(cs1p(:,2)<1,:) = []; cs1p(cs1p(:,2)>nS+1,:) = [];
% 
%     % Subset1 Points to be added are those candidate subset1 points which
%     % are in B but not already in subset1
%     s1p_tba = setdiff( intersect(cs1p, B, 'rows'), ...
%         lcss(1:nr_lcss,:), 'rows');
%     [nr_s1p_tba, ~] = size(s1p_tba);
%     % Put them in subset1
%     lcss(nr_lcss+1:nr_lcss+nr_s1p_tba, :) = s1p_tba;
%     % Update Recently Added Points and the Number of Subset1 Points
%     rap = s1p_tba;
%     nr_lcss = nr_lcss+nr_s1p_tba;
% end
% 
% flag = nr_lcss==nr_B;
% 
% if plot_stuff 
%     figure,
%     plot(y, x, 'kx');
%     axis([0 nS+2 0 nW+2])
%     connectedness_title = ['Level ' num2str(i) ' is '];
%     if flag
%         title([connectedness_title 'connected'], 'FontSize', 20);
%     else
%         title([connectedness_title 'NOT connected'], 'FontSize', 20);
%     end
%     h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
%     h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
%     pause;
%     close
% end
% end
% 
% %% Find largest connected superset
% function lcss = find_largest_connected_superset(nW, nS, B, single_element)
% % Initiate the first subset
% lcss = zeros((nW+1)*(nS+1),2); % Largest Connected Superset
% 
% % Pick the first element of subset1
% lcss(1,:) = single_element;
% nr_lcss = 1;
% rap = lcss(1,:); % Recently Added Points
% while ~isempty(rap)
%     % Find the indices of neighbor points of the recently added points
%     % Candidate Subset1 Points
%     [nr_rap, ~] = size(rap);
%     
%     cs1p = [rap+ones(nr_rap,1)*[-1 0]; rap+ones(nr_rap,1)*[1 0]; ...
%         rap+ones(nr_rap,1)*[0 1]; rap+ones(nr_rap,1)*[0 -1]];
% 
%     % Remove repetitions
%     cs1p = unique(cs1p, 'rows');
%     % Remove out-of-grid points
%     cs1p(cs1p(:,1)<1,:) = []; cs1p(cs1p(:,1)>nW+1,:) = [];
%     cs1p(cs1p(:,2)<1,:) = []; cs1p(cs1p(:,2)>nS+1,:) = [];
%     
%     % Subset1 Points to be added are those candidate subset1 points which
%     % are in B but not already in subset1
%     s1p_tba = setdiff( intersect(cs1p, B, 'rows'), ...
%         lcss(1:nr_lcss,:), 'rows');
%     [nr_s1p_tba, ~] = size(s1p_tba);
%     % Put them in subset1
%     lcss(nr_lcss+1:nr_lcss+nr_s1p_tba, :) = s1p_tba;
%     % Update Recently Added Points and the Number of Subset1 Points
%     rap = s1p_tba;
%     nr_lcss = nr_lcss+nr_s1p_tba;
% end
% 
% % Remove extra zero rows
% lcss(lcss(:,1)==0,:) = [];
% end