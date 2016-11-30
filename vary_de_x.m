function vary_de_x(varargin)
% --------------------------------------------------------------------
plot_specific_density = 1;
% --------------------------------------------------------------------
if strcmp(varargin{1}, 'load from file')
    load(varargin{2});
    fontsize = 25;
    plot_specific_density = 1;
    plot_all_in_one = 1;
    specific_de_idx = 3;
    specific_x_idx = 2;
else
    de = varargin{1};
    x = varargin{2};
    graph_type = varargin{3};
    
    nr_sims = 10;
    nr_events = 1e5;
    % --------------------------------------------------------------------
    exps{length(de),length(x)} = [];
    avg_exp{length(de),length(x)} = [];
    for i=1:length(de)
        for j=1:length(x)
            [exps{i,j}, avg_exp{i,j}] = run_multiple_network_sims(...
                'nr_sims', nr_sims, 'nr_events', nr_events,  ...
                'W', [100 de(i) de(i)/x(j)], 'S', [100 de(i) de(i)/x(j)], ...
                'wake_graph', graph_type, 'sleep_graph', graph_type, ...
                'compute_domains', 0, 'compute_bouts', 0, ...
                'save_workspace', 0);
        end
    end
    
    filename = ['saved_workspaces/various_de_x_' num2str(min(de)) '-' ...
        num2str(max(de)) datestr(now,'_mmmdd_HHMMSS') '.mat'];
    save(filename);
end

% Plot all average density contour plots
if plot_all_in_one
    figure,
    for i=1:length(de)
        for j=1:length(x)
            subplot(length(de), length(x), (i-1)*length(x) + j);
            plot_activity_domains_on_surface_and_contour(avg_exp{i,j}, ...
                'mark_domain_pts', 0, 'plot_on_new_fig', 0, ...
                'want_surface', 0, 'FontSize', fontsize);
            set(gca, 'XTickLabel', '', 'YTickLabel', '', ...
                'xtick', [], 'ytick', []);
            if j==1
                ylabel(['de=' num2str(de(i))], 'FontSize', fontsize);
                %             ylabel(num2str(de(i)/di), 'FontSize', fontsize)
            else
                ylabel('');
            end
            
            if i==1
                title(['x=' num2str(x(j))], 'FontSize', fontsize);
            else
                title('', 'FontSize', fontsize);
            end
            xlabel('');
        end
    end
end

% Plot a specific density on a separate figure
if plot_specific_density
    i = specific_de_idx; %   find(de==specific_de);
    j = specific_x_idx; %      find(di==specific_x);
    %     for j=1:length(di)
    %         for i=1:length(de)
    plot_activity_domains_on_surface_and_contour(avg_exp{i,j}, ...
        'want_contour', 0, 'plot_on_new_fig', 1 , ...
        'mark_domain_pts', 0, 'FontSize', fontsize);
    set(gca, 'XTickLabel', '', 'YTickLabel', '', 'ZTickLabel', '', ...
        'xtick', [], 'ytick', [], 'ztick', []);
    title(avg_exp{i,j}.sim_title(3:17), 'FontSize', fontsize);
    box off;
    %         end
    %     end
end

end