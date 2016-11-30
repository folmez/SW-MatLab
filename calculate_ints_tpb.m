function [duration_matrix, wake_active_domain, sleep_active_domain, ...
    wake_percentage, sleep_percentage] = calculate_ints_tpb(varargin)
% CALCULATE_INTS_TPB calculates the intervals of activity switched between
% two cluster by using a transition path based method. Wake- and
% sleep-active domains are defined as the hottest regions visited by the
% state variable WE and SE. Periods followed by a transition are defined as
% intervals. Duration matrix is the output that contains all intervals
% lengths and types

fprintf('NOTE TO SELF: GO OVER THIS SCRIPT. GOING DOWN\N THE PEAK BY MAKING SURE NEIGHBOR STATES ARE\N GREATER THAN A PERCENTAGE COULD RAISE PROBLEMS')

nW = varargin{1}(1);
nS = varargin{1}(2);
nr_events = varargin{2};
sim_title = varargin{3};
t = varargin{4};
WE = varargin{5};
SE = varargin{6};
A = varargin{7};
display_sim_summary = 1;
switch nargin
    case 8
        display_sim_summary = varargin{8};
end

display_detailed_simulation_summary = 0;
plot_stuff_for_accuracy = 0;
plot_domains_at_every_loop = 0;
plot_pdfs_for_all_intervals = 1;
% ------------------------------------------------------------------------
% KEEP IN MIND THAT A IS (nW+1)x(nS+1) matrix.
% A(1,:)=A(:,1)=0 and A(WE+1,SE+1) is defined as the number of times that
% process visits WE+1,SE+1.

% This function assumes that the matrix A has two peaks!

% When A is plotted using "surf", xlabel is columns, i.e. SE and ylabel is
% rows, i.e. WE

m = 0;
nr_of_means = 3;
% this parameter determines the size scale over which bouts are plotted.
% i.e. tails are not considered!

% quotient determines the # bins used in the bout histograms
quot = 1;

% ------------------------------------------------------------------------
% IDENTIFY THE TWO PEAKS
row_nums = (nW+1:-1:1)'*ones(1,nS+1);
col_nums = ones(nW+1,1)*(1:nS+1);
wake_region = nS*(row_nums-1) < nW*(col_nums-1); % everything above geometric diagonal
sleep_region = nS*(row_nums-1) >= nW*(col_nums-1); % everything below ...
% ... but we need the other diagonal
wake_region = wake_region(:,end:-1:1); % everythin below the other diagonal
sleep_region = sleep_region(:,end:-1:1);
A_wake_part = A; A_wake_part(sleep_region)=0;
A_sleep_part = A; A_sleep_part(wake_region)=0;

% Wake peak
[wpWE,wpSE,~] = find(A_wake_part==max(max(A_wake_part)),1);
wp = [wpWE wpSE];
wpFreq = A(wpWE,wpSE);

% Sleep peak
[spWE,spSE,~] = find(A_sleep_part==max(max(A_sleep_part)),1);
sp = [spWE spSE];
spFreq = A(spWE,spSE);

for delta_WE = -1:1
    for delta_SE = -1:1
        if A(median([1, nW, spWE+delta_WE]), ...
                median([1, nS, spSE+delta_SE])) > spFreq 
            error('Fatih: Sleep domain doesn'' exist!!!');
        elseif A(median([1, nW, wpWE+delta_WE]), ...
                median([1, nS, wpSE+delta_SE])) > wpFreq 
            error('Fatih: Wake domain doesn'' exist!!!');
        end
    end
end

if display_detailed_simulation_summary
    display(wp);
    display(sp);
end
% ------------------------------------------------------------------------
% Find a set of indices of A that are close to wake peak
wake_percentage = 1;
wake_active_domain = wp;
if display_detailed_simulation_summary
    fprintf('WAKE ACTIVE DOMAIN POINTS\nWE+1\tSE+1\n');
end
while size(wake_active_domain,1)==1 || ...
        ~isempty(intersect(wake_active_domain,sp,'rows'))
    wake_percentage = wake_percentage-0.01;
    wake_active_domain = wp;
    nr_good_k_bdry_pts = 1; k=0;
    while nr_good_k_bdry_pts~=0
        nr_good_k_bdry_pts = 0; k = k+1;
        
        % find the indices of the k-th boundary
        bdry_idx = ones(4*(k),1)*wp + [[(-k:k)' ; (k-1:-1:-k+1)' ], ...
            [(0:k-1)';(k:-1:0)' ; (-1:-1:-k)' ; (-k+1:-1)']];
        % remove rows that contain negatives and larger than N+1 indices
        bdry_idx( bdry_idx(:,1)<1 | bdry_idx(:,1)>nW+1 | ...
            bdry_idx(:,2)<1 | bdry_idx(:,2)>nS+1, : ) = [];
        % subscripts of boundary indices to be tested
        good_bdry_idx = sub2ind([nW+1 nS+1],bdry_idx(:,1),bdry_idx(:,2));
        % test boundary points if they occur at least percentage frequently
        good_bdry_idx = good_bdry_idx( A(good_bdry_idx) > ...
            wpFreq*(1-wake_percentage) );
        % record number of good kth boundary points
        nr_good_k_bdry_pts = length(good_bdry_idx);
        % return to 2D indices
        [tempX,tempY] = ind2sub([nW+1 nS+1],good_bdry_idx);
        % recreate a good_bdry_idx
        good_bdry_idx = [tempX,tempY];
        if ~isempty(good_bdry_idx)
            wake_active_domain = union(wake_active_domain, ...
                good_bdry_idx,'rows');
        end
        if display_detailed_simulation_summary
            fprintf('%i\t%i\n',good_bdry_idx');
        end
    end
    wad = wake_active_domain;
    
    if plot_domains_at_every_loop
        figure, contour(A,100); hold,
        plot(wad(:,2),wad(:,1),'kx');
        title([num2str(wake_percentage) ' down the wake peak']); hold off
    end
end
% ------------------------------------------------------------------------
% Find a set of indices of A that are close to wake peak
sleep_percentage = 1;
sleep_active_domain = sp;
while size(sleep_active_domain,1)==1 || ...
        ~isempty(intersect(sleep_active_domain,wp,'rows'))
    sleep_percentage = sleep_percentage-0.01;
    sleep_active_domain = sp;
    nr_good_k_bdry_pts = 1; k=0;
    while nr_good_k_bdry_pts~=0
        nr_good_k_bdry_pts = 0; k = k+1;
        
        % find the indices of the k-th boundary
        bdry_idx = ones(4*(k),1)*sp + [[(-k:k)' ; (k-1:-1:-k+1)' ], ...
            [(0:k-1)';(k:-1:0)' ; (-1:-1:-k)' ; (-k+1:-1)']];
        % remove rows that contain negatives and larger than N+1 indices
        bdry_idx( bdry_idx(:,1)<1 | bdry_idx(:,1)>nW+1 | ...
            bdry_idx(:,2)<1 | bdry_idx(:,2)>nS+1, : ) = [];
        % subscripts of boundary indices to be tested
        good_bdry_idx = sub2ind([nW+1 nS+1],bdry_idx(:,1),bdry_idx(:,2));
        % test boudnary points if they occur at least percentage frequently
        good_bdry_idx = good_bdry_idx( A(good_bdry_idx) > ...
            spFreq*(1-sleep_percentage));
        % record number of good kth boundary points
        nr_good_k_bdry_pts = length(good_bdry_idx);
        % return to 2D indices
        [tempX,tempY] = ind2sub([nW+1 nS+1],good_bdry_idx);
        % recreate a good_bdry_idx
        good_bdry_idx = [tempX,tempY];
        if ~isempty(good_bdry_idx)
            sleep_active_domain = union(sleep_active_domain, ...
                good_bdry_idx,'rows');
        end
    end
    sad = sleep_active_domain;
    
    if plot_domains_at_every_loop
        figure, contour(A,100); hold,
        plot(sad(:,2),sad(:,1),'kx');
        title([num2str(sleep_percentage) ' down the sleep peak']); hold off
    end
end
% ------------------------------------------------------------------------
% Plot wake-active and sleep-active domains on surface
figure, surf(A), hold,
plot3(wad(:,2),wad(:,1),A(sub2ind([nW+1 nS+1],wad(:,1),wad(:,2))),'kx', ...
    wpSE,wpWE,A(sub2ind([nW+1 nS+1],wpWE,wpSE)),'bx');
plot3(sad(:,2),sad(:,1),A(sub2ind([nW+1 nS+1],sad(:,1),sad(:,2))),'ko', ...
    spSE,spWE,A(sub2ind([nW+1 nS+1],spWE,spSE)),'bo');
xlabel('SE'); ylabel('WE');title(sim_title);
legend('PDF of SE x WE','Wake-active domain','Sleep-active domain', ...
    'Location','Best');


% Plot wake-active and sleep-active domains on contour
figure, contour(A,50); hold,
plot(wad(:,2),wad(:,1),'kx');
plot(sad(:,2),sad(:,1),'ko');
xlabel('SE'); ylabel('WE'); title(sim_title);
legend('PDF of SE x WE','Wake-active domain','Sleep-active domain', ...
    'Location','Best');
% ------------------------------------------------------------------------
if ~isempty(intersect(wad,sad,'rows'))
    error('Wake and sleep domains overlap!!!');
end
% ------------------------------------------------------------------------

newState = zeros(nr_events,1);
% Note that the matrix A is defined by adding +1 to WE and SE!
tAD = tic;
for i=1:length(wad)
    newState = newState + (WE==wad(i,1)-1 & SE==wad(i,2)-1);
end

for i=1:length(sad)
    newState = newState - (WE==sad(i,1)-1 & SE==sad(i,2)-1);
end
fprintf('\nDetermining activity domains took %3.2f minutes\n',toc(tAD)/60);

% Compute wake, sleep, non-transition and real-transition bouts
durations = diff([0;t([find(diff(newState)~=0);nr_events])]);
duration_types = newState([find(diff(newState)~=0);nr_events]);

% If first and last durations happen to be "transitions", reassign them
if duration_types(1)==0, duration_types(1)=duration_types(2); end
if duration_types(end)==0, duration_types(end)=duration_types(end-1); end

% Replace every 0 for transitions, with 0.5*(previous bout + next bouts).
% This makes sure that 1 0 1 --> 1 1 1 and -1 0 -1 --> -1 -1 -1 and 1 0 -1,
% -1 0 1 --> * 0 *
idx = find(duration_types==0);
duration_types(idx) = 0.5*(duration_types(idx-1)+duration_types(idx+1));

% Now recalculate durations with the updated wake and sleep durations
new_t = cumsum(durations);
newState = duration_types;
new_nr_events = length(duration_types);
new_durations = diff([0;new_t([find(diff(newState)~=0);new_nr_events])]);
new_duration_types = newState([find(diff(newState)~=0);new_nr_events]);

new_t = cumsum(new_durations);
newState = new_duration_types;

si = new_durations(new_duration_types==-1);
wi = new_durations(new_duration_types==1);
ti = new_durations(new_duration_types==0);

duration_matrix = [new_durations new_duration_types];
wake_active_domain = wad-1;
sleep_active_domain = sad-1;

% ------------------------------------------------------------------------
bt = ti(ti<nr_of_means*mean(ti) & ti>m);
[btn,btx] = hist(bt,round(max(bt)-min(bt))/quot);

figure,
subplot(1,2,1); loglog(btx,btn,'b.-');
title('Log-log and Semi-log plots of trans bouts')
subplot(1,2,2); semilogy(btx,btn,'b.-');
title(sim_title);
% ------------------------------------------------------------------------
if plot_stuff_for_accuracy
    for i=0:10
        figure
        hold
        
        X1 = i*1e4+1:(i+1)*1e4;
        plot(t(X1),WE(X1),'b',t(X1),SE(X1),'g');
        
        X2 = find(new_t < t((i+1)*1e4) & new_t > t(i*1e4+1));
        for j=1:length(X2)
            plot([new_t(X2(j)) new_t(X2(j))],[5 95],'r');
        end
        
        hold off
    end
end
% ------------------------------------------------------------------------
if plot_pdfs_for_all_intervals
    figure, subplot(1,3,1);
    plot_approximate_pdf_of_data(si, 'data_title', 'Sleep intervals', ...
        'plot_on_new_fig', 0);
    subplot(1,3,2);
    plot_approximate_pdf_of_data(wi, 'data_title', 'Wake intervals', ...
        'plot_on_new_fig', 0);
    subplot(1,3,3);
    plot_approximate_pdf_of_data(ti, 'data_title', 'Transition intervals', ...
        'plot_on_new_fig', 0);
end
% ------------------------------------------------------------------------
if display_sim_summary
    fprintf('\nINTERVAL CALCULATION SUMMARY\n');
    fprintf(sim_title);
    fprintf('\nNumber of sleep intervals: %i',length(si));
    fprintf('\nNumber of wake intervals: %i',length(wi));
    fprintf('\nNumber of transition intervals: %i',length(ti));
    
    fprintf('\nMean sleep interval: %6.4f',mean(si));
    fprintf('\nMean wake interval: %6.4f',mean(wi));
    fprintf('\nMean transition interval: %6.4f',mean(ti));
    
    fprintf('\nFraction of sleep time: %6.4f',sum(si)/t(end));
    fprintf('\nFraction of wake time: %6.4f',sum(wi)/t(end));
    fprintf('\nFraction of transition time: %6.4f\n',sum(ti)/t(end));
end
end