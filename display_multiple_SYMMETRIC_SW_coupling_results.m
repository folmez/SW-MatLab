function display_multiple_SYMMETRIC_SW_coupling_results(varargin)
avg_exp = varargin{1};
%-------------------------------------------------------------------------
fprintf('\n\n\tCOMBINED SLEEP WAKE BOUTS for symmetric graph choices (Half transitions)\n\n');
for k = 1:length(avg_exp)
    wake = avg_exp{k}.bouts.wake.half_trans;
    sleep = avg_exp{k}.bouts.sleep.half_trans;
    display_combined_bout_results
end
%-------------------------------------------------------------------------
fprintf('\n\n\tCOMBINED SLEEP WAKE BOUTS for symmetric graph choices (No transitions)\n\n');
for k = 1:length(avg_exp)
    wake = avg_exp{k}.bouts.wake.no_trans;
    sleep = avg_exp{k}.bouts.sleep.no_trans;
    display_combined_bout_results
end
%-------------------------------------------------------------------------
    function display_combined_bout_results
        if k==1
            fprintf('Sleep graph\t\tWake graph\t\t');
            fprintf('alpha (std)\tvalid(%%)\t');
            fprintf('beta\n');
        end
        temp_graph_titles = avg_exp{k}.graph_titles;
        all_alphas = [wake.PL.alpha_vector; sleep.PL.alpha_vector];
        all_betas = [wake.exp_tail.beta_vector; sleep.exp_tail.beta_vector];
        fprintf('%s\t%s\t', temp_graph_titles.sleep, temp_graph_titles.wake);
        fprintf('%1.2f (%1.2f)\t%2.0i\t\t', ...
            mean(all_alphas), std(all_alphas), ...
            100*length(all_alphas)/(2*avg_exp{k}.nr_sims));
        fprintf('%1.1e\n', mean(all_betas));
    end
%-------------------------------------------------------------------------
end