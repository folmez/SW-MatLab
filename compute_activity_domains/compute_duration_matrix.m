function duration_matrix = compute_duration_matrix(varargin)
addpath ../../power-law-estimator/

%   COMPUTE_DURATION_MATRIX calculates intervals of activity switched
%   between two graphs by using a transition-based method. Duration matrix
%   is the output that contains all interval lengths and types

%% Input parameters
t = varargin{1};
WE = varargin{2};
SE = varargin{3};
wad = varargin{4};
sad = varargin{5};
nr_events = varargin{6};
sim_title = varargin{7};
display_duration_matrix_summary = 1;
plot_transition_bouts = 0;
plot_stuff_for_accuracy = 0;
plot_pdfs_for_all_intervals = 0;

i = 8;
while i<=length(varargin),
    switch varargin{i},
        case 'plot_pdfs_for_all_intervals'
            plot_pdfs_for_all_intervals = varargin{i+1};
        case 'plot_stuff_for_accuracy'
            plot_stuff_for_accuracy = varargin{i+1};
        case 'plot_transition_bouts'
            plot_transition_bouts = varargin{i+1};
        case 'display_duration_matrix_summary'
            display_duration_matrix_summary = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

%% Model
newState = zeros(nr_events,1);
% Note that the matrix A is defined by adding +1 to WE and SE!
tAD = tic;
[lwad, ~] = size(wad);
for i=1:lwad
    newState = newState + (WE==wad(i,1) & SE==wad(i,2));
end
[lsad, ~] = size(sad);
for i=1:lsad
    newState = newState - (WE==sad(i,1) & SE==sad(i,2));
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

si = new_durations(new_duration_types==-1);
wi = new_durations(new_duration_types==1);
ti = new_durations(new_duration_types==0);

duration_matrix = [new_durations new_duration_types];

%% Plot stuff
if plot_transition_bouts
    bt = ti(ti<nr_of_means*mean(ti) & ti>m);
    [btn,btx] = hist(bt,round(max(bt)-min(bt))/quot);
    figure,
    subplot(1,2,1); loglog(btx,btn,'b.-');
    title('Log-log and Semi-log plots of trans bouts')
    subplot(1,2,2); semilogy(btx,btn,'b.-');
    title(sim_title);
end

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

if plot_pdfs_for_all_intervals
    figure,
    subplot(1,2+~isempty(ti),1);
    plot_approximate_pdf_of_data(si, 'data_title', 'Sleep intervals', ...
        'plot_on_new_fig', 0);
    subplot(1,2+~isempty(ti),2);
    plot_approximate_pdf_of_data(wi, 'data_title', 'Wake intervals', ...
        'plot_on_new_fig', 0);
    if ~isempty(ti)
        subplot(1,3,3);
        plot_approximate_pdf_of_data(ti, 'data_title', ...
            'Transition intervals', 'plot_on_new_fig', 0);
    end
end

%% Display stuff
if display_duration_matrix_summary
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