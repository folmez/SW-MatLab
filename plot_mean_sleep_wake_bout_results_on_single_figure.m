function plot_mean_sleep_wake_bout_results_on_single_figure(varargin)
exps = varargin{1};
% ------------------------------------------------------------------------
mean_bout_matrix = zeros(length(exps), 6);  % W - no, prev, half
                                            % S - no, prev, half
for i=1:length(exps)
    mean_bout_matrix(i,1) = mean(exps(i).bouts.wake.no_trans.data);
    mean_bout_matrix(i,2) = mean(exps(i).bouts.wake.prev_trans.data);    
    mean_bout_matrix(i,3)= mean(exps(i).bouts.wake.half_trans.data);    
    mean_bout_matrix(i,4)= mean(exps(i).bouts.sleep.no_trans.data);    
    mean_bout_matrix(i,5)= mean(exps(i).bouts.sleep.prev_trans.data);    
    mean_bout_matrix(i,6)= mean(exps(i).bouts.sleep.half_trans.data);
end
% ------------------------------------------------------------------------
% Display mean bout results
fprintf('MEAN\tW-no\tW-prev\tW-half\tS-no\tS-prev\tS-half\n');
fprintf('MEAN\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f\n', mean_bout_matrix');
fprintf('\n');
fprintf('MEAN\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f\t%3.1f\n', mean(mean_bout_matrix));
% ------------------------------------------------------------------------
end