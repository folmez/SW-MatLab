function avg_exp = display_sleep_wake_bouts_exponential_tail_results(varargin)
exps = varargin{1};
display_individual_exponential_tail_results = varargin{2};
avg_exp = varargin{3};
% ------------------------------------------------------------------------
% No transition bouts
no_trans = zeros(length(exps), 5);  % Data size,
                                    % W-beta, xmin
                                    % S-beta, xmin
for i=1:length(exps)
    no_trans(i,:) = [length(exps(i).bouts.wake.no_trans.data), ...
        exps(i).bouts.wake.no_trans.exp_tail.beta, ...
        exps(i).bouts.wake.no_trans.exp_tail.xmin, ...
        exps(i).bouts.sleep.no_trans.exp_tail.beta, ...
        exps(i).bouts.sleep.no_trans.exp_tail.xmin];
end

% Previous transitions are added
prev_trans = zeros(length(exps), 5);
for i=1:length(exps)
    prev_trans(i,:) = [length(exps(i).bouts.wake.prev_trans.data), ...
        exps(i).bouts.wake.prev_trans.exp_tail.beta, ...
        exps(i).bouts.wake.prev_trans.exp_tail.xmin, ...
        exps(i).bouts.sleep.prev_trans.exp_tail.beta, ...
        exps(i).bouts.sleep.prev_trans.exp_tail.xmin];
end

% Half transition bouts
half_trans = zeros(length(exps), 5);
for i=1:length(exps)
    half_trans(i,:) = [length(exps(i).bouts.wake.half_trans.data), ...
        exps(i).bouts.wake.half_trans.exp_tail.beta, ...
        exps(i).bouts.wake.half_trans.exp_tail.xmin, ...
        exps(i).bouts.sleep.half_trans.exp_tail.beta, ...
        exps(i).bouts.sleep.half_trans.exp_tail.xmin];
end
% ------------------------------------------------------------------------
if display_individual_exponential_tail_results
    fprintf('EXPONENTIAL TAILS\n');
    fprintf('\t\t\t\t\tWAKE BOUTS\t\t\t\tSLEEP BOUTS\n\n');
    
    fprintf('NO TRANS');
    display_bout_exp_tail_matrix(no_trans);
    
    fprintf('\nPREV TRANS');
    display_bout_exp_tail_matrix(prev_trans);
    
    fprintf('\nHALF TRANS');
    display_bout_exp_tail_matrix(half_trans);
end
% ---
    function display_bout_exp_tail_matrix(bout_exp_tails)
        fprintf('\tSize');
        fprintf('\t\tbeta\txmin');
        fprintf('\t\tbeta\txmin\n');
        for k = 1:length(exps)
            fprintf('\t\t%1.1e', bout_exp_tails(k, 1));
            fprintf('\t\t%1.1e\t%3.2f', bout_exp_tails(k, 2:3));
            fprintf('\t\t%1.1e\t%3.2f', bout_exp_tails(k, 4:5));
            fprintf('\n');
        end
    end
% ------------------------------------------------------------------------
fprintf('\nAVERAGES of ');
fprintf([exps(1).sim_title '\n\n']);
fprintf('EXPONENTIAL TAILS\n');
fprintf('\t\t\t\tWAKE BOUTS\t\tSLEEP BOUTS\n\n');

fprintf('NO TRANS');
[avg_exp.bouts.wake.no_trans.exp_tail, avg_exp.bouts.sleep.no_trans.exp_tail] = ...
    display_MEAN_STD_bout_exp_tail_matrix(no_trans);

fprintf('PREV TRANS');
[avg_exp.bouts.wake.prev_trans.exp_tail, avg_exp.bouts.sleep.prev_trans.exp_tail] = ...
    display_MEAN_STD_bout_exp_tail_matrix(prev_trans);

fprintf('HALF TRANS');
[avg_exp.bouts.wake.half_trans.exp_tail, avg_exp.bouts.sleep.half_trans.exp_tail] = ...
    display_MEAN_STD_bout_exp_tail_matrix(half_trans);
% ---
    function [wake_exp_tail, sleep_exp_tail] = ...
            display_MEAN_STD_bout_exp_tail_matrix(bout_PLs)
        
        fprintf('\tSize');
        fprintf('\t\tbeta\txmin');
        fprintf('\t\tbeta\txmin\n');
        
        fprintf('MEAN\t\t%1.1e', mean(bout_PLs(:, 1),1));
        fprintf('\t\t%1.1e\t%3.2f', mean(bout_PLs(:, 2:3),1));
        fprintf('\t\t%1.1e\t%3.2f', mean(bout_PLs(:, 4:5),1));
        fprintf('\n');
        
        fprintf('STD\t\t%1.1e', std(bout_PLs(:, 1),1));
        fprintf('\t\t%1.1e\t%3.2f', std(bout_PLs(:, 2:3),1));
        fprintf('\t\t%1.1e\t%3.2f', std(bout_PLs(:, 4:5),1));
        fprintf('\n\n');
        
        wake_exp_tail_avg = mean(bout_PLs(:, 2:3),1);
        sleep_exp_tail_avg = mean(bout_PLs(:, 4:5),1);
        
        wake_exp_tail.beta = wake_exp_tail_avg(1);
        wake_exp_tail.xmin = wake_exp_tail_avg(2);
        wake_exp_tail.beta_vector = bout_PLs(:, 2);
        
        sleep_exp_tail.beta = sleep_exp_tail_avg(1);
        sleep_exp_tail.xmin = sleep_exp_tail_avg(2);
        sleep_exp_tail.beta_vector = bout_PLs(:, 4);
   end
% ------------------------------------------------------------------------
% avg_exp.bouts.wake.no_trans.exp_tail.beta = wake_nt_exp_tail_avg(1);
% avg_exp.bouts.wake.no_trans.exp_tail.xmin = wake_nt_exp_tail_avg(2);
% avg_exp.bouts.sleep.no_trans.exp_tail.beta = sleep_nt_exp_tail_avg(1);
% avg_exp.bouts.sleep.no_trans.exp_tail.xmin = sleep_nt_exp_tail_avg(2);
% 
% avg_exp.bouts.wake.prev_trans.exp_tail.beta = wake_pt_exp_tail_avg(1);
% avg_exp.bouts.wake.prev_trans.exp_tail.xmin = wake_pt_exp_tail_avg(2);
% avg_exp.bouts.sleep.prev_trans.exp_tail.beta = sleep_pt_exp_tail_avg(1);
% avg_exp.bouts.sleep.prev_trans.exp_tail.xmin = sleep_pt_exp_tail_avg(2);
% 
% 
% avg_exp.bouts.wake.half_trans.exp_tail.beta = wake_ht_exp_tail_avg(1);
% avg_exp.bouts.wake.half_trans.exp_tail.xmin = wake_ht_exp_tail_avg(2);
% avg_exp.bouts.sleep.half_trans.exp_tail.beta = sleep_ht_exp_tail_avg(1);
% avg_exp.bouts.sleep.half_trans.exp_tail.xmin = sleep_ht_exp_tail_avg(2);
end