function display_multiple_SW_coupling_results(varargin)
avg_exps = varargin{1};
%-------------------------------------------------------------------------
fprintf('\t\t\t\tWAKE BOUTS (Half transitions)\n\n');
fprintf('Sleep graph\t\tWake graph\t\talpha (std)\txmin\txmax\tvalid\t\tbeta\txmin\n');
for k = 1:length(avg_exps)
    temp_graph_titles = avg_exps{k}.graph_titles;
    temp_PL = avg_exps{k}.bouts.wake.half_trans.PL;
    temp_exp_tail = avg_exps{k}.bouts.wake.half_trans.exp_tail;
    fprintf('%s\t%s\t%1.2f (%1.2f)\t%3.2f\t%3.2f\t%1.2f\t\t%1.1e\t%3.2f\n', ...
        temp_graph_titles.sleep, temp_graph_titles.wake, ...
        temp_PL.alpha, temp_PL.alpha_std, temp_PL.xmin, temp_PL.xmax, ...
        temp_PL.acceptance_ratio, ...
        temp_exp_tail.beta, temp_exp_tail.xmin);
end

fprintf('\n\t\t\t\tWAKE BOUTS (No transitions)\n\n');
fprintf('Sleep graph\t\tWake graph\t\talpha (std)\txmin\txmax\tvalid\t\tbeta\txmin\n');
for k = 1:length(avg_exps)
    temp_graph_titles = avg_exps{k}.graph_titles;
    temp_PL = avg_exps{k}.bouts.wake.no_trans.PL;
    temp_exp_tail = avg_exps{k}.bouts.wake.no_trans.exp_tail;
    fprintf('%s\t%s\t%1.2f (%1.2f)\t%3.2f\t%3.2f\t%1.2f\t\t%1.1e\t%3.2f\n', ...
        temp_graph_titles.sleep, temp_graph_titles.wake, ...
        temp_PL.alpha, temp_PL.alpha_std, temp_PL.xmin, temp_PL.xmax, ...
        temp_PL.acceptance_ratio, ...
        temp_exp_tail.beta, temp_exp_tail.xmin);
end

fprintf('\n\t\t\t\tWAKE BOUTS (Previous transitions)\n\n');
fprintf('Sleep graph\t\tWake graph\t\talpha (std)\txmin\txmax\tvalid\t\tbeta\txmin\n');
for k = 1:length(avg_exps)
    temp_graph_titles = avg_exps{k}.graph_titles;
    temp_PL = avg_exps{k}.bouts.wake.prev_trans.PL;
    temp_exp_tail = avg_exps{k}.bouts.wake.prev_trans.exp_tail;
    fprintf('%s\t%s\t%1.2f (%1.2f)\t%3.2f\t%3.2f\t%1.2f\t\t%1.1e\t%3.2f\n', ...
        temp_graph_titles.sleep, temp_graph_titles.wake, ...
        temp_PL.alpha, temp_PL.alpha_std, temp_PL.xmin, temp_PL.xmax, ...
        temp_PL.acceptance_ratio, ...
        temp_exp_tail.beta, temp_exp_tail.xmin);
end

fprintf('\n\t\t\t\tSLEEP BOUTS (Half transitions)\n\n');
fprintf('Sleep graph\t\tWake graph\t\talpha (std)\txmin\txmax\tvalid\t\tbeta\txmin\n');
for k = 1:length(avg_exps)
    temp_graph_titles = avg_exps{k}.graph_titles;
    temp_PL = avg_exps{k}.bouts.sleep.half_trans.PL;
    temp_exp_tail = avg_exps{k}.bouts.sleep.half_trans.exp_tail;
    fprintf('%s\t%s\t%1.2f (%1.2f)\t%3.2f\t%3.2f\t%1.2f\t\t%1.1e\t%3.2f\n', ...
        temp_graph_titles.sleep, temp_graph_titles.wake, ...
        temp_PL.alpha, temp_PL.alpha_std, temp_PL.xmin, temp_PL.xmax, ...
        temp_PL.acceptance_ratio, ...
        temp_exp_tail.beta, temp_exp_tail.xmin);
end

fprintf('\n\t\t\t\tSLEEP BOUTS (No transitions)\n\n');
fprintf('Sleep graph\t\tWake graph\t\talpha (std)\txmin\txmax\tvalid\t\tbeta\txmin\n');
for k = 1:length(avg_exps)
    temp_graph_titles = avg_exps{k}.graph_titles;
    temp_PL = avg_exps{k}.bouts.sleep.no_trans.PL;
    temp_exp_tail = avg_exps{k}.bouts.sleep.no_trans.exp_tail;
    fprintf('%s\t%s\t%1.2f (%1.2f)\t%3.2f\t%3.2f\t%1.2f\t\t%1.1e\t%3.2f\n', ...
        temp_graph_titles.sleep, temp_graph_titles.wake, ...
        temp_PL.alpha, temp_PL.alpha_std, temp_PL.xmin, temp_PL.xmax, ...
        temp_PL.acceptance_ratio, ...
        temp_exp_tail.beta, temp_exp_tail.xmin);
end

fprintf('\n\t\t\t\tSLEEP BOUTS (Previous transitions)\n\n');
fprintf('Sleep graph\t\tWake graph\t\talpha (std)\txmin\txmax\tvalid\t\tbeta\txmin\n');
for k = 1:length(avg_exps)
    temp_graph_titles = avg_exps{k}.graph_titles;
    temp_PL = avg_exps{k}.bouts.sleep.prev_trans.PL;
    temp_exp_tail = avg_exps{k}.bouts.sleep.prev_trans.exp_tail;
    fprintf('%s\t%s\t%1.2f (%1.2f)\t%3.2f\t%3.2f\t%1.2f\t\t%1.1e\t%3.2f\n', ...
        temp_graph_titles.sleep, temp_graph_titles.wake, ...
        temp_PL.alpha, temp_PL.alpha_std, temp_PL.xmin, temp_PL.xmax, ...
        temp_PL.acceptance_ratio, ...
        temp_exp_tail.beta, temp_exp_tail.xmin);
end
%-------------------------------------------------------------------------
fprintf('\t\t\t\tWAKE BOUTS (Half transitions)\n\n');
fprintf('Sleep graph\t\tWake graph\t\t#events\tMean #bouts\tMean bout size\n');
for k = 1:length(avg_exps)
    temp_graph_titles = avg_exps{k}.graph_titles;
    temp_PL = avg_exps{k}.bouts.wake.half_trans.PL;
    fprintf('%s\t%s\t%1.1e\t%1.1e\t\t%3.2f\n', ...
        temp_graph_titles.sleep, temp_graph_titles.wake, ...
        avg_exps{k}.nr_events, avg_exps{k}.mean_nr_bouts, temp_PL.mean_bout_size);
end
%-------------------------------------------------------------------------
end