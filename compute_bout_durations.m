function varargout = compute_bout_durations(varargin)
%% Input arguments
duration_matrix = varargin{1};
compute_bouts = varargin{2};

%% Model
if compute_bouts
    [sleep_bouts_prev_trans, wake_bouts_prev_trans, ~] = ...
        SW_interval_to_bout(duration_matrix, ...
        'transition intervals are added to previous bout');
    [sleep_intervals, wake_intervals, transition_intervals] = ...
        SW_interval_to_bout(duration_matrix, 'intervals are bouts');
    [sleep_bouts_half_trans, wake_bouts_half_trans, ~] = ...
        SW_interval_to_bout(duration_matrix, ...
        ['transition intervals are cut in half and ' ...
        'added to closest wake/sleep interval']);
else
    sleep_bouts_prev_trans = [];
    wake_bouts_prev_trans = [];
    sleep_intervals = [];
    wake_intervals = [];
    transition_intervals = [];
    sleep_bouts_half_trans = [];
    wake_bouts_half_trans = [];
end

%% Output arguments
varargout{1} = sleep_bouts_prev_trans;
varargout{2} = wake_bouts_prev_trans;
varargout{3} = sleep_intervals;
varargout{4} = wake_intervals;
varargout{5} = transition_intervals;
varargout{6} = sleep_bouts_half_trans;
varargout{7} = wake_bouts_half_trans;

end


%% Transform intervals to bouts based on different definitions

function [sb, wb, tb, wstb, swtb] = ...
    SW_interval_to_bout(duration_matrix, which_tpb_method)
% Input:    Duration matrix (first column is sizes of bouts, second column
%                               is types of bouts)
%           Which transition method (add or not add transitions:
%                                       'intervals are bouts'
%                                       'transition intervals are added to
%                                        previous bout'
%                                       'transition intervals are cut in 
%                                        Half and added to closest 
%                                        wake/sleep interval'
%
% Output:       sb: sleep bouts
%               wb: wake bouts
%               tb: transition bouts

display_sim_summary = 1;

durs = duration_matrix(:,1); % durations
types = duration_matrix(:,2); % duration types
total_time = sum(durs);
% Note that calculate_ints_tpb generates duration_matrix so that first and
% last types are not transitions. Therefore length of types is an odd
% number
% ------------------------------------------------------------------------
% which_tpb_method could be a string or a number
if ~ischar(which_tpb_method)
    if which_tpb_method == 1
        which_tpb_method = 'intervals are bouts';
    elseif which_tpb_method == 2
        which_tpb_method = ['transition intervals are added ' ...
            'to previous bout'];
    elseif which_tpb_method==3
        which_tpb_method = ['transition intervals are cut in ' ...
            'half and added to closest wake/sleep interval'];
    else
        error('Wrong input in turning intervals to bouts!');
    end
end

% Calculate bouts
if strcmp(which_tpb_method, 'intervals are bouts')
    sb = durs(types==-1);
    wb = durs(types==1);
    tb = durs(types==0);
    
    tb_idx = find(types==0);
    prev_bouts = types(tb_idx-1);
    
    % wstb: wake to sleep transition bouts
    % swtb: sleep to wake transition bouts
    wstb_idx = tb_idx(prev_bouts==1);
    wstb = durs(wstb_idx);
    swtb_idx = tb_idx(prev_bouts==-1);
    swtb = durs(swtb_idx);
elseif strcmp(which_tpb_method, ...
        'transition intervals are added to previous bout')
    idx = find(types==0);
    durs(idx-1) = durs(idx-1)+durs(idx);
    durs(idx) = [];
    types(idx) = [];
    sb = durs(types==-1);
    wb = durs(types==1);
    tb = []; wstb = []; swtb = [];
elseif strcmp(which_tpb_method, ['transition intervals are cut in ' ...
        'half and added to closest wake/sleep interval'])
    idx = find(types==0);
    durs(idx-1) = durs(idx-1)+0.5*durs(idx);
    durs(idx+1) = durs(idx+1)+0.5*durs(idx);
    durs(idx) = [];
    types(idx) = [];
    sb = durs(types==-1);
    wb = durs(types==1);
    tb = []; wstb = []; swtb = [];
end

% Display simulation summary
if display_sim_summary
    fprintf('\nINTERVAL CALCULATION SUMMARY\n');
    fprintf(which_tpb_method);
    fprintf('\nNumber of sleep bouts: %i',length(sb));
    fprintf('\nNumber of wake bouts: %i',length(wb));
    fprintf('\nNumber of transition bouts: %i', length(tb));
    
    fprintf('\nMean sleep bout: %6.4f', mean(sb));
    fprintf('\nMean wake bout: %6.4f', mean(wb));
    fprintf('\nMean transition bout: %6.4f', mean(tb));
    fprintf('\nMean awake->asleep transition bout: %6.4f', mean(wstb));
    fprintf('\nMean asleep->awake transition bout: %6.4f', mean(swtb));
    
    fprintf('\nFraction of sleep time: %6.4f', sum(sb)/total_time);
    fprintf('\nFraction of wake time: %6.4f', sum(wb)/total_time);
    fprintf('\nFraction of transition time: %6.4f\n', sum(tb)/total_time);
end

end