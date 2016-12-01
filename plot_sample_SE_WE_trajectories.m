function plot_sample_SE_WE_trajectories(varargin)
%% Input arguments
duration_matrix = varargin{1};
wake_active_domain = varargin{2};
sleep_active_domain = varargin{3};
SE = varargin{4};
WE = varargin{5};
sim_title = varargin{6};
nS = varargin{7};
nW = varargin{8};
plot_sample_trajectories = varargin{9};

%% Plot sample trajectories
if plot_sample_trajectories
    make_movie_from_bout_trajs = 0;
    what_type_of_bouts = [-1];
    
    min_bout_size = 0;
    max_bout_size = inf;
    nr_bout_trajs = 50;
    pause_after_figure = 1;
    
    if make_movie_from_bout_trajs
        close all;
    end
    
    for i=1:nr_bout_trajs
        bout_index_set = find(duration_matrix(:,1) > min_bout_size & ...
            duration_matrix(:,1) < max_bout_size);
        bout_index = bout_index_set(randperm(length(bout_index_set),1));
        while ~ismember(duration_matrix(bout_index,2),what_type_of_bouts)
            bout_index = bout_index_set( ...
                randperm(length(bout_index_set),1));
        end
        
        figure, hold on;
        plot(wake_active_domain(:,2), wake_active_domain(:,1), 'bo');
        plot(sleep_active_domain(:,2), sleep_active_domain(:,1), 'gd');
        
        for j = 0:2
            bout_begin_time = sum(duration_matrix(1:bout_index+j-1,1));
            bout_end_time = bout_begin_time + ...
                duration_matrix(bout_index+j,1);
            bout_type = duration_matrix(bout_index+j,2);
            
            traj_index_set_min = find(t >= bout_begin_time,1);
            traj_index_set_max = find(t > bout_end_time,1);
            traj_index_set = (traj_index_set_min:traj_index_set_max-1)';
            if bout_type==-1
                str_bout_type = 'Sleep';
                plot(SE(traj_index_set(1)), WE(traj_index_set(1)), 'g*');
                plot(SE(traj_index_set), WE(traj_index_set), 'g.-');
            elseif bout_type==1
                str_bout_type = 'Wake';
                plot(SE(traj_index_set(1)), WE(traj_index_set(1)), 'b*');
                plot(SE(traj_index_set), WE(traj_index_set), 'b.-');
            elseif bout_type==0
                str_bout_type = 'Transition';
                plot(SE(traj_index_set(1)), WE(traj_index_set(1)), 'k*');
                plot(SE(traj_index_set), WE(traj_index_set), 'k.-');
            end
            %plot(SE(traj_index_set(1)), WE(traj_index_set(1)), 'r*');
            %plot(SE(traj_index_set(end)), WE(traj_index_set(end)), 'k*');
        end
        legend(sim_title,'Location','Best');
        axis([0 nS+1 0 nW+1]); xlabel('SE'); ylabel('WE');
        title(['Sample ' str_bout_type ' bout, size = ' ...
            num2str(duration_matrix(bout_index,1))]);
        
        if pause_after_figure, pause; end
    end
end

end