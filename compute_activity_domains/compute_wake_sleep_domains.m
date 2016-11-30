function [wad, sad] = compute_wake_sleep_domains(varargin)
%   [wad, sad] = COMPUTE_WAKE_SLEEP_DOMAINS([nW nS], sim_title, A);
%   COMPUTE_WAKE_SLEEP_DOMAINS computes wake and sleep domains for a given
%   simulated (WE(t),SE(t)). An error is displayed if there is only one
%   state.
% ------------------------------------------------------------------------
nW = varargin{1}(1);
nS = varargin{1}(2);
sim_title = varargin{2};
A = varargin{3};
% ------------------------------------------------------------------------
i = 8;
while i<=length(varargin),
    switch varargin{i},
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
make_quick_existence_of_two_domain_check = 1;
enter_peaks_manually = 0;
plot_stable_eq_points = 0;
% ------------------------------------------------------------------------
% IDENTIFY WAKE and SLEEP PEAKS
% Find split state space (grid) into two triangle-like regions and define
% peaks in each region accordingly as wake or sleep peak.
row_nums = (nW+1:-1:1)'*ones(1,nS+1);
col_nums = ones(nW+1,1)*(1:nS+1);
wake_reg = nS*(row_nums-1) < nW*(col_nums-1); % everything above geometric diagonal
sleep_reg = nS*(row_nums-1) >= nW*(col_nums-1); % everything below ...
% ... but we need the other diagonal
wake_reg = wake_reg(:,end:-1:1); % everything below the other diagonal
sleep_reg = sleep_reg(:,end:-1:1);
A_wake_part = A; A_wake_part(sleep_reg)=0;
A_sleep_part = A; A_sleep_part(wake_reg)=0;

if enter_peaks_manually
    % Wake peak
    wp = input('Enter the wake peak as [SE,WE]: ');
    % Sleep peak
    sp = input('Enter the sleep peak as [SE,WE]: ');    
else
    % Wake peak
    [wp(1), wp(2), ~] = find(A_wake_part==max(max(A_wake_part)),1);
    % Sleep peak
    [sp(1), sp(2), ~] = find(A_sleep_part==max(max(A_sleep_part)),1);
end

% Make a quick check whether everthing is OK
if make_quick_existence_of_two_domain_check
    for delta_WE = -1:1
        for delta_SE = -1:1
            if A(median([1, nW, sp(1)+delta_WE]), ...
                    median([1, nS, sp(2)+delta_SE])) > A(sp(1),sp(2))
                warning('Fatih: Sleep domain doesn'' exist!!!');
            elseif A(median([1, nW, wp(1)+delta_WE]), ...
                    median([1, nS, wp(2)+delta_SE])) > A(wp(1),wp(2))
                warning('Fatih: Wake domain doesn'' exist!!!');
            end
        end
    end
end
% ------------------------------------------------------------------------
% Determine wake and sleep domains as the largest connected supersets
% containing one peak and not containing the other

wad = find_largest_domain(nW, nS, A, wp, sp);
sad = find_largest_domain(nW, nS, A, sp, wp);

% Plot wake- and sleep- active domain
figure, contour(A,50); hold,
plot(wad(:,2),wad(:,1),'kx');
plot(sad(:,2),sad(:,1),'ko');
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

% Plot the two stable equilibrium points that can be calculated for the
% two-state model - INCOMPLETE
if plot_stable_eq_points
    dW = varargin{4}(1);
    dS = varargin{4}(2);
    dIW = varargin{4}(3);
    dIS = varargin{4}(4);
    lambda_e = varargin{4}(5);
    lambda_i = varargin{4}(6);
    
    x = mean([dW dS])/mean([dIW dIS]);
    y = lambda_e/lambda_i;
    
    C = 1/2-(x+1)./(2*x*(y-1));
    D = sqrt((1-x).*(x*y-1).*(x*y-x.^2*y-3*x-1)) ./ (2*x*(y-1).*(1-x));
end
end
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
function ad = find_largest_domain(nW, nS, A, peak1, peak2)
%   This function determines an active domain as the largest connected
%   superset containing peak1 and not containing peak2
ad = peak1;
A_vals = unique(sort(A(:)));
A_vals(A_vals>A(peak1(1),peak1(2))) = [];
i = length(A_vals);
while 1
    % Find next larger superlevel set
    [x,y,~] = find(A>=A_vals(i));
    % Make sure next superlevel set does contain the peak point #1
    if ~ismember(peak1, [x y], 'rows')
        error('A superlevel set didn''t contain peak point!');
    end
    % Determine largest connected superset of peak point #1
    ad_next = find_largest_connected_superset(nW, nS, [x y], peak1);
    % If the above set contains peak point #2, break the while loop and
    % return previous active domain. Otherwise, go on to the next larger
    % superlevel set
    if ismember(peak2, ad_next, 'rows')
        break;
    else
        ad = ad_next;
        i = i-1;
    end
end
end
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
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