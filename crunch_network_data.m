function crunch_network_data(varargin) %which_filename, which_exps)
addpath ../test_rbm_data/
addpath ../power-law-estimator/

% what_to_do = 'run power-law fits and p-values';
% what_to_do = 'bout_plots_for_DS13_poster';
% what_to_do = 'alpha, beta errorbars for DS13 poster';
% what_to_do = 'Plot power-law fits on bout distributions';
% what_to_do = 'Plot mean wake, sleep and transition times';
what_to_do = 'Plot average densities';

if strcmp(what_to_do,'run power-law fits and p-values')
    %% Run Power-law fits and p-values to network data
    filename{1} = 'SF_7pnt6_6pnt8_SF_7pnt6_6pnt8_five_sims.mat';
    filename{2} = 'ER_7pnt6_6pnt8_ER_7pnt6_6pnt8_three_sims.mat';
    filename{3} = 'ER_100_3pnt8_4_ER_100_3pnt8_4_five_sims.mat';
    filename{4} = 'SF_100_3pnt8_4pnt5_SF_100_3pnt8_4pnt5_five_sims.mat';
    filename{5} = 'ER_15_20_ER_15_20_one_sim.mat';
    filename{6} = 'SF_100_3pnt8_4pnt5_SF_100_3pnt8_4pnt5_four_sims_round2.mat';
    filename{7} = 'ER_100_3pnt9_4pnt5_SF_100_3pnt9_4pnt5_five_sims.mat';
    filename{8} = 'ER_50_ER_300_four_sims.mat';
    filename{9} = 'ER_20_10_5_ER_100_2pnt6_0pnt85_six_sims.mat';
    
    which_filename = varargin{1}; which_exps = varargin{2};
    
    crunch_transition_bouts = 0;
    crunch_wake_bouts = 0;
    crunch_sleep_bouts = 0;
    save_figures = 0;
    need_p_value = 0;
    plot_stuff = 0;
    
    load(['saved_workspaces/' filename{which_filename}]);
    
    for jj = which_exps
        experiment =  exps(jj);
        
        plot_activity_domains_on_surface_and_contour('exp',experiment);
        
        int_to_bout_defn{1} = 'intervals are bouts';
        int_to_bout_defn{2} = 'transition intervals are added to previous bout';
        int_to_bout_defn{3} = 'transition intervals are cut in half and added to closest wake/sleep interval';
        
        sim_title = experiment.sim_title;
        display(sim_title);
        
        for k = 1
            [Sb, Wb, Tb] = interval_to_bout(experiment.duration_matrix, ...
                int_to_bout_defn{k});
            
            if length(Wb)>1e5
                warning('Number of bouts is too high! Consider only first 50k!!!');
                Wb = Wb(1:1e5); Sb = Sb(1:1e5); Tb = Tb(1:1e5);
            end
            
            plot_approximate_pdf_of_data(Wb(1:1e3), 'data_title', 'Wake bouts');
            pause;
            find_pl_fit_with_adapt_slope_detect(Wb(1:1e3), ...
                'data_title', sim_title,'nr_reps',5);
            
            % WAKE BOUTS
            if crunch_wake_bouts
                [wake_pl, wake_exp_tail] = find_exp_tail_and_pl_fit( ...
                    Wb, sim_title, ...
                    'data_title', ...
                    ['Wake bouts' num2str(length(Wb),' (#bouts=%2.1e') ...
                    num2str(mean(Wb),', mean wake bout size=%4.2f)')], ...
                    'nr_trial_points', nr_trial_points, ...
                    'need_p_value' , need_p_value, ...
                    'qof_methods', [1, 2, 5], ...
                    'which_pl_fit_on_plot', 5, ...
                    'plot_stuff', plot_stuff, ...
                    'KS_slope', 0.002, ... % 02, ...
                    'search_pl_in_tailless_data', 1, ...
                    'search_pl_beyond_this', 0);
                
                % pl_or_exp_better_fit(Wb, wake_pl.xmin, wake_pl.xmax);
            end
            
            % Approximate PDF of the transition bouts
            figure, 
            if crunch_transition_bouts && k==1
                [Tbn, Tbx] = hist(Tb,logspace(log10(min(Tb)/2),log10(max(Tb)+1),1000));
                Tbn = Tbn./diff([0 Tbx]);
                loglog(Tbx,Tbn);
                xlabel(['Transition bouts' ...
                    num2str(mean(Tb),'(mean transition bout size=%4.2f)')]);
                title(sim_title);
            end
            
            % SLEEP BOUTS
            if crunch_sleep_bouts
                [sleep_pl, sleep_exp_tail] = find_exp_tail_and_pl_fit( ...
                    Sb, sim_title, ...
                    'data_title', ...
                    ['Sleep bouts' num2str(length(Sb),' (%2.1e'), ...
                    num2str(mean(Sb),', mean sleep bout size=%4.2f)')], ...
                    'nr_trial_points', nr_trial_points, 'need_p_value' , 0, ...
                    'qof_methods', [1, 2, 5], ...
                    'which_pl_fit_on_plot', 5, 'plot_stuff', plot_stuff, ...
                    'KS_slope', 0.002, ... % 02, ...
                    'search_pl_in_tailless_data', 1, ...
                    'search_pl_beyond_this', 0);
                
                % pl_or_exp_better_fit(Sb, sleep_pl.xmin, sleep_pl.xmax);
            end
            
            
            if save_figures
                mkdir('figures'); cd figures;
                mkdir(filename{which_filename}(1:end-4)); % just filename without .mat
                cd(filename{which_filename}(1:end-4))
                mkdir(['sim' num2str(jj)]); cd(['sim' num2str(jj)])
                figure(1);
                saveas(gcf, ['WE_SE_surf_sim' num2str(jj)], 'fig');
                saveas(gcf, ['WE_SE_surf_sim' num2str(jj)], 'png');
                
                figure(2);
                saveas(gcf, ['WE_SE_contour_sim' num2str(jj)], 'fig');
                saveas(gcf, ['WE_SE_contour_sim' num2str(jj)], 'png');
                
                figure(3);
                saveas(gcf, ['wake_bouts_sim' num2str(jj)], 'fig');
                saveas(gcf, ['wake_bouts_sim' num2str(jj)], 'png');
                
                figure(4);
                saveas(gcf, ['trans_ints_sim' num2str(jj)], 'fig');
                saveas(gcf, ['trans_ints_sim' num2str(jj)], 'png');
                
                figure(5);
                saveas(gcf, ['sleep_bouts_sim' num2str(jj)], 'fig');
                saveas(gcf, ['sleep_bouts_sim' num2str(jj)], 'png');
                
                close all; cd .., cd .., cd ..
            end
            
        end
    end
    
    % BELOW IS NEEDED TO WRAP UP WORKSPACES
%         load('two_cluster_ER_20_10_5.05_ER_100_2.64_0.78_May11_001242.mat');
%     exps(6) = experiment; clearvars -except exps
%     load('two_cluster_ER_20_10_5.05_ER_100_2.64_0.82_May11_001638.mat');
%     exps(5) = experiment; clearvars -except exps
%         load('two_cluster_ER_20_10_5.05_ER_100_2.64_0.82_May11_001947.mat');
%         exps(4) = experiment; clearvars -except exps
%         load('two_cluster_ER_20_10_5_ER_100_2.56_0.83_May11_001412.mat');
%         exps(3) = experiment; clearvars -except exps
%         load('two_cluster_ER_20_10_5_ER_100_2.62_0.76_May11_002246.mat');
%         exps(2) = experiment; clearvars -except exps
%         load('two_cluster_ER_20_10_5_ER_100_2.6_0.85_May11_002126.mat');
%         exps(1) = experiment; clearvars -except exps
%         save('ER_20_10_5_ER_100_2pnt6_0pnt85_six_sims.mat','exps');

elseif strcmp(what_to_do,'Plot power-law fits on bout distributions')
    %%
    clear all
    filename = 'two_cluster_Mar29_0840_p_values.mat';
    load(filename);
    
    % Plot power-law fits
    figure,
    % Plot sleep bounds and the sleep power-law fit
    subplot(2,1,1), loglog(Sbx,Sbn); title(sim_title);
    xlabel('Sleep bout sizes'); ylabel('Frequency');
    sleep_pl_lower_cutoff = sleep_pl_bounds(5,2);
    sleep_pl_upper_cutoff = sleep_pl_bounds(5,3);
    sleep_pl_exponent = (-1)*sleep_pl_bounds(5,1);
    sleep_pl_KS = sleep_pl_bounds(5,5);
    
    hold on;
    x0 = sleep_pl_lower_cutoff; y0 = Sbn(find(Sbx>x0,1));
    x1 = sleep_pl_upper_cutoff; slope = sleep_pl_exponent;
    plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'ro-','Linewidth',2);
    legend('Approximate PDF of sleep bouts',...
        ['Fitted (using two-step) power-law in [' num2str(x0) ',' num2str(x1) ...
        '] w/ alpha=' num2str(slope)], ...
        'Location','Best');
    hold off;
    
    % Plot wake bounds and the wake power-law fit
    subplot(2,1,2), loglog(Wbx,Wbn);
    title(sim_title);
    xlabel('Wake bout sizes'); ylabel('Frequency');
    wake_pl_lower_cutoff = wake_pl_bounds(5,2);
    wake_pl_upper_cutoff = wake_pl_bounds(5,3);
    wake_pl_exponent = (-1)*wake_pl_bounds(5,1);
    wake_pl_KS = wake_pl_bounds(5,5);
    
    hold on;
    x0 = wake_pl_lower_cutoff; y0 = Wbn(find(Wbx>x0,1));
    x1 = wake_pl_upper_cutoff; slope = wake_pl_exponent;
    plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'ro-','Linewidth',2);
    hold off;
    legend('Approximate PDF of wake bouts',...
        ['Fitted (using two-step) power-law in [' num2str(x0) ',' num2str(x1) ...
        '] w/ alpha=' num2str(slope)], ...
        'Location','Best');
    
    % Plot transition intervals
    [Tbn, Tbx] = hist(Tb,round((max(Tb)-min(Tb))/quotient));
    figure,
    semilogy(Tbx, Tbn); xlabel('Transition intervals'); title(sim_title);
    
elseif strcmp(what_to_do,'bout_plots_for_DS13_poster')
    %% Plot certain bouts and power-law fits on the same figure
    nr_trial_points = 100;
    save_figures = 0;
    plot_stuff = 0;

    which_filename = (1:4)';
    which_experiment_and_bout = [3 1; 1 1; 1 -1; 2 -1]; % 1 is Wake, -1 is Sleep bout
    
    for jj = 1:length(which_filename)
        load(filename{which_filename(jj)});
        experiment =  exps(which_experiment_and_bout(jj,1));
        [Sb, Wb, ~] = interval_to_bout(experiment.duration_matrix, ...
            'intervals are bouts');
        plot_activity_domains_on_surface_and_contour('exp',experiment);
    
        sim_title = experiment.sim_title;
        display(sim_title);

        if which_experiment_and_bout(jj,2) == 1 % WAKE BOUTS
            [wake_pl, wake_exp_tail] = find_exp_tail_and_pl_fit( ...
                Wb, sim_title, ...
                'data_title', ...
                ['Wake bouts' num2str(length(Wb),' (#bouts=%2.1e') ...
                num2str(mean(Wb),', mean wake bout size=%4.2f)')], ...
                'nr_trial_points', nr_trial_points, 'need_p_value' , 0, ...
                'qof_methods', [1, 2, 5], ...
                'which_pl_fit_on_plot', 5, 'plot_stuff', plot_stuff, ...
                'KS_slope', 0.002, ... % 02, ...
                'search_pl_in_tailless_data', 1, ...
                'search_pl_beyond_this', 0);
            X{jj} = Wb; pl(jj) = wake_pl; exp_tail(jj) = wake_exp_tail;
        elseif which_experiment_and_bout(jj,2) == -1 % SLEEP BOUTS
            [sleep_pl, sleep_exp_tail] = find_exp_tail_and_pl_fit( ...
                Sb, sim_title, ...
                'data_title', ...
                ['Sleep bouts' num2str(length(Sb),' (%2.1e'), ...
                num2str(mean(Sb),', mean sleep bout size=%4.2f)')], ...
                'nr_trial_points', nr_trial_points, 'need_p_value' , 0, ...
                'qof_methods', [1, 2, 5], ...
                'which_pl_fit_on_plot', 5, 'plot_stuff', plot_stuff, ...
                'KS_slope', 0.002, ... % 02, ...
                'search_pl_in_tailless_data', 1, ...
                'search_pl_beyond_this', 0);
            X{jj} = Sb; pl(jj) = sleep_pl; exp_tail(jj) = sleep_exp_tail;
        end
        simulation_title{jj} = sim_title;
    end
    
    cc = jet(length(which_filename)); % Color map
    
    % Plot all power-laws on one figure
    figure, hold on
    for jj = 1:4%[1, 2, 4, 3]
        binning = 'log';
        curve_shift = 10^(1-jj);
        data_set = X{jj}(X{jj}<exp_tail(jj).xmin);
        if strcmp(binning,'log')
            [Xn, Xx] = hist(data_set, ...
                logspace(log10(min(data_set)/2),log10(max(data_set)+1),500));
            Xn = Xn./diff([0 Xx]);
        elseif strcmp(binning,'linear')
            quotient = 0.1;
            [Xn, Xx] = hist(data_set,10000); % round((max(X)-min(X))/quotient));
        end
        
        % Approximate PDF of data
        plot(Xx*curve_shift,Xn*curve_shift, 'color', cc(jj,:));
        
        % Power-law fit
        x0 = pl(jj).xmin*curve_shift;
        y0 = Xn(find(Xx*curve_shift>x0,1))*curve_shift;
        x1 = pl(jj).xmax*curve_shift; slope = pl(jj).alpha;
        plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'color', cc(jj,:),'Linewidth',2);
    end
    legend(simulation_title{1}, pl(1).title, ...
        simulation_title{2}, pl(2).title, ...
        simulation_title{3}, pl(3).title, ...
        simulation_title{4}, pl(4).title, 'Location', 'Best')
    set(gca,'xscale','log'); set(gca,'yscale','log');
    set(gca,'xticklabel',[]); set(gca,'yticklabel',[])
    xlabel('Bout size'); ylabel('Bout distribution');
    axis tight;
    hold off;
    
    
    % Plot all exponential tails on one figure
    figure, hold on;
    for jj=1:length(which_filename)
        binning = 'log';
        curve_shift = 10^(1-jj);
%         if jj==1, curve_shift = 1e-5; elseif jj==2, curve_shift=1e-2; end
        data_set = X{jj}(X{jj}>pl(jj).xmax);
        if strcmp(binning,'log')
            [Xn, Xx] = hist(data_set, ...
                logspace(log10(min(data_set)/2),log10(max(data_set)+1),400));
            Xn = Xn./diff([0 Xx]);
        elseif strcmp(binning,'linear')
            quotient = 0.1;
            [Xn, Xx] = hist(data_set,500); % round((max(X)-min(X))/quotient));
        end
        
        % Approximate PDF of data
        plot(Xx,Xn*curve_shift, 'color', cc(jj,:));
        
        x0 = exp_tail(jj).xmin;
        y0 = Xn(find(Xx>x0,1))*curve_shift;
        slope = exp_tail(jj).beta;
        %         plot(sort(data_set(data_set>=x0)), ...
        %             y0*exp(slope*(sort(data_set(data_set>=x0))-x0)), 'Linewidth',2);
        plot([x0 max(data_set)/2],...
            [y0,y0*exp(slope*(max(data_set)/2-x0))], ...
            'color', cc(jj,:), 'Linewidth',2);
    end
    legend(simulation_title{1}, exp_tail(1).title, ...
        simulation_title{2}, exp_tail(2).title, ...
        simulation_title{3}, exp_tail(3).title, ...
        simulation_title{4}, exp_tail(4).title, 'Location', 'Best')
    set(gca,'yscale','log');
    set(gca,'xticklabel',[]); set(gca,'yticklabel',[])
    xlabel('Bout size'); ylabel('Bout distribution');
    axis tight;
    hold off;
elseif strcmp(what_to_do, 'alpha, beta errorbars for DS13 poster')
    %% alpha and beta that are on the pictures
    ER_ER_low = [0.98 0.981 5.2e-4 2.6e-3 ; 0.887 1.032 4.7e-4 6.4e-4; ...
        0.937 0.89 6.1e-4 5.5e-4; 1.025 1.024 1.9e-3 2.6e-4];
    SF_SF_low = [0.946 0.886 1.5e-3 2e-3; 0.998 0.826 1.4e-3 9.5e-4; ...
        0.647 0.922 5.4e-4 1.8e-3; 1.165 0.737 7.1e-4 1.2e-4];
    ER_ER_high = [0.66 0.668 6.1e-4 4.8e-4; 0.704 0.678 5e-4 5.8e-4; ...
        0.642 0.647 6.4e-4 6.7e-4];
    SF_SF_high = [0.65 0.568 2.1e-3 1.8e-3; 0.61 0.655 2.5e-3 1.2e-3; ...
        0.876 0.582 2.5e-3 1.2e-3; 0.596 0.612 1.4e-3 1.4e-3];
    
    data_set{4} = SF_SF_low;
    data_set{3} = ER_ER_low;
    data_set{2} = ER_ER_high;
    data_set{1} = SF_SF_high;
    
    cc = jet(4); % Color map
    
    % Power-law exponents on errorbars
    figure, subplot(2,1,1);
    hold on;
    for jj=1:4
        alphas = reshape(data_set{jj}(:,1:2), ...
            numel(data_set{jj}(:,1:2)), 1);
        errorbar(jj, mean(alphas), std(alphas), ...
            'Color', cc(jj,:), 'Linewidth', 2);
    end
    h_legend = legend('SF(100,7.6,6.8)', 'ER(100,7.6,6.8)', ...
        'ER(100,3.8,4.5)', 'SF(100.3.8,4.5)','Location','Best');
    set(h_legend,'FontSize',13);
    set(gca,'xticklabel',[]); ylabel('alphas');
    h_ylabel = get(gca,'YLabel'); set(h_ylabel,'FontSize',15);
    title('Power-law exponents for multiple random realizations', ...
        'FontSize',15);
    grid on;
    hold off;
    
    % Exponent decay rates on errorbars
    subplot(2,1,2);
    hold on;
    for jj=1:4
        betas = reshape(data_set{jj}(:,3:4), ...
            numel(data_set{jj}(:,3:4)), 1);
        errorbar(jj, mean(betas), std(betas), ...
            'Color', cc(jj,:), 'Linewidth', 2);
    end
%     h_legend = legend('SF(100,7.6,6.8)', 'ER(100,7.6,6.8)', ...
%         'ER(100,3.8,4.5)', 'SF(100.3.8,4.5)','Location','Best');
%     set(h_legend,'FontSize',13);
    set(gca,'yscale','log');
    set(gca,'xticklabel',[]); ylabel('betas');
    h_ylabel = get(gca,'YLabel'); set(h_ylabel,'FontSize',15);
    title('Exponential tail decay rates for multiple random realizations', ...
        'FontSize',15);
    grid on;
    hold off
    
elseif strcmp(what_to_do,'Plot mean wake, sleep and transition times')
    %% Plot all mean times for each network coupling

    nr_filenames = length(filename);
    figure,
    set(gca,'yscale','log');
    title('Mean sleep, wake, transition intervals', 'FontSize', 15);
    hold on;
    for i = 1:nr_filenames
        load(['saved_workspaces/' filename{i}]);
        
        nr_exps = length(exps);
        mean_SWT = zeros(nr_exps, 3); % mean sleep, wake, transition bouts
        std_SWT = zeros(nr_exps, 3); % std sleep, wake, transition bouts
        sim_title{i} = exps(1).sim_title;
        for j = 1:nr_exps
            experiment = exps(j);

            % plot_activity_domains_on_surface_and_contour('exp',experiment);
            
            [Sb, Wb, Tb] = interval_to_bout(experiment.duration_matrix, ...
                'intervals are bouts');
            
            mean_SWT(j,:) = [mean(Sb) mean(Wb) mean(Tb)];
            std_SWT(j,:) = [std(Sb) std(Wb) std(Tb)];
        end

        semilogy(ones(nr_exps,1) + 0.1*i, mean_SWT(:,1), 'g*-', ...
            1 + 0.1*i, mean(mean_SWT(:,1)), 'go');
        semilogy(2*ones(nr_exps,1) + 0.1*i, mean_SWT(:,2),'b*-', ...
            2 + 0.1*i, mean(mean_SWT(:,2)), 'bo');
        semilogy(3*ones(nr_exps,1) + 0.1*i, mean_SWT(:,3),'k*-', ...
            3 + 0.1*i, mean(mean_SWT(:,3)), 'ko');
        
        text(4.5, 8e3/2^i, [num2str(i) ') ' sim_title{i}]);
    end
    xlim([0.9 6.5]);
elseif strcmp(what_to_do,'Plot average densities')
    %% Plot average densities to show stochastic bifurcation from one state to two states
    
%     filename{1} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_3_4.5_ER_100_3_4.5_Jan22_171705.mat';
%     filename{2} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_4_4.5_ER_100_4_4.5_Jan22_170636.mat';
%     filename{3} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_5_4.5_ER_100_5_4.5_Jan22_170542.mat';
%     filename{4} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_6_4.5_ER_100_6_4.5_Jan22_170611.mat';
%     filename{5} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_7_4.5_ER_100_7_4.5_Jan22_170512.mat';
%     filename{6} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_8_4.5_ER_100_8_4.5_Jan22_170538.mat';
    
    filename{1} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_3_4.5_ER_100_3_4.5_Jan22_191933.mat';
    filename{2} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_4_4.5_ER_100_4_4.5_Jan22_191505.mat';
    filename{3} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_5_4.5_ER_100_5_4.5_Jan22_190539.mat';
    filename{4} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_6_4.5_ER_100_6_4.5_Jan22_185234.mat';
    filename{5} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_7_4.5_ER_100_7_4.5_Jan22_184314.mat';
    filename{6} = 'saved_workspaces/average_density_plots/mult_network_100_sims_ER_100_8_4.5_ER_100_8_4.5_Jan22_183849.mat';
    
    
    figure,
    for i=length(filename):-1:1
        avg_exp(i) = run_multiple_network_sims(...
            'load_from_file', filename{i}, ...
            'plot_avg_density', 0);
        
        subplot(3,2,i);
        plot_activity_domains_on_surface_and_contour(avg_exp(i), ...
        'mark_domain_pts', 0, 'want_contour', 0, 'plot_on_new_fig', 0);
    end
end