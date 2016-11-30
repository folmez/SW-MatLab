function avg_exp = display_sleep_wake_bouts_PL_results(varargin)
exps = varargin{1};
display_individual_PL_results = varargin{2};
avg_exp = varargin{3};
% ------------------------------------------------------------------------
% No transition bouts
no_trans = zeros(length(exps), 11); % Data size,
                                    % W-alpha, xmin, xmax, xmax/xmin
                                    % S-alpha, xmin, xmax, xmax/xmin
                                    % Mean wake bout size,
                                    % Mean sleep bout size
for i=1:length(exps)
    wake = exps(i).bouts.wake.no_trans;
    sleep = exps(i).bouts.sleep.no_trans;
    no_trans(i,:) = put_all_wake_sleep_data_in_a_row(wake, sleep);
end

% Previous transitions are added
prev_trans = zeros(length(exps), 11);
for i=1:length(exps)
    wake = exps(i).bouts.wake.prev_trans;
    sleep = exps(i).bouts.sleep.prev_trans;
    prev_trans(i,:) = put_all_wake_sleep_data_in_a_row(wake, sleep);
end

% Half transition bouts
half_trans = zeros(length(exps), 11);
for i=1:length(exps)
    wake = exps(i).bouts.wake.half_trans;
    sleep = exps(i).bouts.sleep.half_trans;
    half_trans(i,:) = put_all_wake_sleep_data_in_a_row(wake, sleep);
end
% ------------------------------------------------------------------------
    function out_row = put_all_wake_sleep_data_in_a_row(wake, sleep)
        out_row = [length(wake.data), ...
            wake.PL.alpha, ...
            wake.PL.xmin, ...
            wake.PL.xmax, ...
            wake.PL.xmax/wake.PL.xmin, ...
            sleep.PL.alpha, ...
            sleep.PL.xmin, ...
            sleep.PL.xmax, ...
            sleep.PL.xmax/sleep.PL.xmin, ...
            mean(wake.data), ...
            mean(sleep.data)];
    end
% ------------------------------------------------------------------------
if display_individual_PL_results
    fprintf('\t\t\t\t\tWAKE BOUTS\t\t\t\tSLEEP BOUTS\n\n');
    
    fprintf('NO TRANS');
    display_bout_PL_matrix(no_trans);
    
    fprintf('\nPREV TRANS');
    display_bout_PL_matrix(prev_trans);
    
    fprintf('\nHALF TRANS');
    display_bout_PL_matrix(half_trans);
end
% ------------------------------------------------------------------------
    function display_bout_PL_matrix(bout_PLs)
        fprintf('\tSize');
        fprintf('\t\talpha\txmin\txmax\tPL-lgth');
        fprintf('\t\talpha\txmin\txmax\tPL-lgth\n');
        [nr_sims, ~] = size(bout_PLs);
        for k = 1:nr_sims
            fprintf('\t\t%1.1e', bout_PLs(k, 1));
            fprintf('\t\t%1.2f\t%3.2f\t%3.2f\t%1.1e', bout_PLs(k, 2:5));
            fprintf('\t\t%1.2f\t%3.2f\t%3.2f\t%1.1e', bout_PLs(k, 6:9));
            fprintf('\n');
        end
    end
% ------------------------------------------------------------------------
fprintf('\nAVERAGES of ');
fprintf([exps(1).sim_title '\n\n']);
fprintf('POWER-LAWs\n');
fprintf('\t\t\t\t\tWAKE BOUTS\t\t\t\tSLEEP BOUTS\n\n');

fprintf('NO TRANS');
[avg_exp.bouts.wake.no_trans.PL, avg_exp.bouts.sleep.no_trans.PL] = ...
    display_MEAN_STD_bout_PL_matrix(no_trans);

fprintf('PREV TRANS');
[avg_exp.bouts.wake.prev_trans.PL, avg_exp.bouts.sleep.prev_trans.PL] = ...
    display_MEAN_STD_bout_PL_matrix(prev_trans);

fprintf('HALF TRANS');
[avg_exp.bouts.wake.half_trans.PL, avg_exp.bouts.sleep.half_trans.PL] = ...
    display_MEAN_STD_bout_PL_matrix(half_trans);
% ------------------------------------------------------------------------
    function [wake_bout_PL, sleep_bout_PL] = ...
            display_MEAN_STD_bout_PL_matrix(bout_PLs)
        idx1 = ~~bout_PLs(:,2);
        idx2 = ~~bout_PLs(:,6);
        
        fprintf('\tSize');
        fprintf('\t\talpha\txmin\txmax\tPL-lgth');
        fprintf('\t\talpha\txmin\txmax\tPL-lgth\n');
        
        fprintf('MEAN\t\t%1.1e', mean(bout_PLs(idx1, 1),1));
        fprintf('\t\t%1.2f\t%3.2f\t%3.2f\t%1.1e', ...
            mean(bout_PLs(idx1, 2:5),1));
        fprintf('\t\t%1.2f\t%3.2f\t%3.2f\t%1.1e', ...
            mean(bout_PLs(idx2, 6:9),1));
        fprintf('\n');
        
        fprintf('STD\t\t%1.1e', std(bout_PLs(idx1, 1),1));
        fprintf('\t\t%1.2f\t%3.2f\t%3.2f\t%1.1e', ...
            std(bout_PLs(idx1, 2:5),1));
        fprintf('\t\t%1.2f\t%3.2f\t%3.2f\t%1.1e', ...
            std(bout_PLs(idx2, 6:9),1));
        fprintf('\n\n');
        
        wake_PL = [mean(bout_PLs(idx1, 2:4),1) std(bout_PLs(idx1, 2),1) ...
            sum(idx1)/length(bout_PLs(:,2)) mean(bout_PLs(:,10))];
        sleep_PL = [mean(bout_PLs(idx2, 6:8),1) std(bout_PLs(idx2, 6),1) ...
            sum(idx2)/length(bout_PLs(:,6)) mean(bout_PLs(:,11))];
        
        wake_bout_PL.alpha = wake_PL(1);
        wake_bout_PL.xmin = wake_PL(2);
        wake_bout_PL.xmax = wake_PL(3);
        wake_bout_PL.alpha_std = wake_PL(4);
        wake_bout_PL.acceptance_ratio = wake_PL(5);
        wake_bout_PL.mean_bout_size = wake_PL(6);
        wake_bout_PL.alpha_vector = bout_PLs(idx1, 2);
        
        sleep_bout_PL.alpha = sleep_PL(1);
        sleep_bout_PL.xmin = sleep_PL(2);
        sleep_bout_PL.xmax = sleep_PL(3);
        sleep_bout_PL.alpha_std = sleep_PL(4);
        sleep_bout_PL.acceptance_ratio = sleep_PL(5);
        sleep_bout_PL.mean_bout_size = sleep_PL(6);
        sleep_bout_PL.alpha_vector = bout_PLs(idx2, 6);
    end
% ------------------------------------------------------------------------
avg_exp.nr_events = exps(1).nr_events;
avg_exp.nr_sims = length(exps);
avg_exp.mean_nr_bouts = mean(no_trans(:, 1));
end