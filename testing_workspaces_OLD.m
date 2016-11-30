% exp1:
% Wake cluster = SF(100,9.51), Sleep cluster = SF(100,9.51)
% mean number of inhbitory links = 63.63
% number of events = 1e7

% exp2:
% Wake cluster = Erdos-Renyi(100,10.26), Sleep cluster = Erdos-Renyi(100,9.88)
% mean number of inhbitory links = 62.21
% number of events = 1e7

% exp3:
% Wake cluster = Erdos-Renyi(100,4.4), Sleep cluster = Erdos-Renyi(100,4)
% mean number of inhbitory links = 4.01
% number of events = 6e7

% exp4:
% Wake cluster = SF(100,3.85), Sleep cluster = SF(100,3.86)
% mean number of inhbitory links = 4.27
% number of events = 6e7

% SECOND ROUND (making sure of bistability in the heat maps and asymmetricity)

% exp1:
% filename = two_cluster_Jul07_0810.mat
% Wake cluster = ER(100,3,4.5), Sleep cluster = ER(100,4.2,3.3)
% number of events = 1e8
% wake_domain_radius_percnt = 0.60;
% sleep_domain_radius_percnt = 0.51;
%

% exp2:
% filename = two_cluster_Jul27_1721.mat
% Wake cluster = ER(100,5,5), Sleep cluster = SF(100,5.77,4.25)
% number of events = 1e8
% wake_domain_radius_percnt = 0.39;
% sleep_domain_radius_percnt = 0.57;
%

% exp3:
% filename = two_cluster_Jul07_1819.mat
% Wake cluster = ER(100,10,63), Sleep cluster = ER(100,15,45)
% number of events = 1e8
% wake_domain_radius_percnt = 0.76;
% sleep_domain_radius_percnt = 0.85;
%

% exp4:
% Wake cluster = ER(100,5,5), Sleep cluster = ER(100,5.77,4.25)
% number of events = 1e7
% wake_domain_radius_percnt = 0.67;
% sleep_domain_radius_percnt = 0.80;
%

% exp5:
% Wake cluster = ER(100,8,5.8), Sleep cluster = ER(100,7.2,6.5)
% number of events = 1e8
% wake_domain_radius_percnt = 0.99;
% sleep_domain_radius_percnt = 0.99;
%

% % THIRD ROUND (1e7)
% W=ER(100,3,6), S=ER(100,4,5), # of events=10000000
%
% W=ER(100,3,6), S=ER(100,6,3.25), # of events=10000000
%
% W=ER(100,40,40), S=ER(100,30,50), # of events=10000000

% Fourth Round
% calculate_ints_tpb('two_cluster_Aug13_2043.mat',0.89,0.98)
% calculate_ints_tpb('two_cluster_Aug14_0850.mat',0.82,0.96)
% calculate_ints_tpb('two_cluster_Aug14_0449.mat',0.94,0.97)
% calculate_ints_tpb('two_cluster_Aug14_0047.mat',0.99,0.99)

% ----------------------------------------------------------

% averaging_distinct_networks_same_params.mat
% ER-ER && ER-SF && SF-SF experiments (8 each = 24 total)

% ----------------------------------------------------------

% m = 1000;
%
% nr_of_means = 3000000;
%
% this parameter determines the size scale over which bouts are plotted.
% i.e. tails are not considered!

% quotient determines the # bins used in the bout histograms
% EXP# structures are in bouts.mat
%---------------------------------------------------------------

which_bout_definition = 'transition path based';
% which_bout_definition = 'standard';
% 'standard' or 'threshold' or 'ignore noisy bouts' or 'transition path based'
% or 'compare two definition' (then i=j)

which_tpb_method = 'intervals are bouts';
% 'intervals are bouts'
% 'transition intervals are added to previous bout'
% 'transition intervals are cut in half and added to closest wake/sleep interval'

wake_domain_radius_percnt = 0.39;
sleep_domain_radius_percnt = 0.57;
% needed when "transition path based"
% any trajectory that enters neighborhood of a peak that occurs at leat %40
% likely as the peak is considered to be assigned to the state of the peak

max_noisy_bout_length = 5; % needed when "ignore noisy bouts"
% bouts that last less than period are added to the previous bout

threshold = 0.5;
% needed when "threshold"
% awake period = From WE > 0.7*N to SE > 0.7*N
% asleep period = From SE > 0.7*N to WE > 0.7*N

% OPTION 1
plot_individual_experiment_results = 0;
use_current_workspace = 0;

% OPTION 2
compare_all_of_current_round_experiments = 0;

% OPTION 3
plot_sample_trajectories_during_bouts = 0;

% OPTION 4
compare_two_experiments = 0;
% compare experiments i and j
i=1; j=2;

% OPTION 5
compare_tpb_methods_and_standard = 0;

% OPTION 6
plot_drift_diffusion_peclet = 0;

% OPTION 7
compare_transition_intervals = 0;

% OPTION 8
compare_statistics = 0;

% OPTION 9
compute_SW_WS_transition_probability_densities = 0;

% OPTION 10
plot_WE_SE_trajectories_vs_t = 0;

% OPTION 11
plot_sleep_wake_survivals_and_WS_SW_transition_cdfs = 0;

% OPTION 12
adaptive_power_law_detection = 0;

% OPTION 13
detect_power_law_exponent_by_varying_min_bout_size = 0;

% OPTION 14
detect_power_law_exponent_automatically = 1;

plot_sample_simulations = 0;
X = 1:1e4; % Number of events that we want to show

if plot_individual_experiment_results
    which_experiment = 4;
    plot_full = 1;
    plot_heat_map = 0;
    
    use_clauset = 0;
    use_linear_regression = 1;
    
    % min_bout = 1;
    % nr_of_means = 1;
    max_bout = inf;
    quotient1 = 10;
    quotient2 = 100;
    
    bout_definition = which_bout_definition;
    tic
    if ~use_current_workspace
        sleep_bouts = experiments(which_experiment).sleep_bouts;
        wake_bouts = experiments(which_experiment).wake_bouts;
        
        sim_title = experiments(which_experiment).info;
        A = experiments(which_experiment).A;
        WE = experiments(which_experiment).WE;
        SE = experiments(which_experiment).SE;
        nr_events = experiments(which_experiment).nr_events;
        t = experiments(which_experiment).t;
        
        [tempSb,tempWb,tempTb] = interval_to_bout(...
            experiments(which_experiment).duration_matrix,which_tpb_method);
        
        % if plot_heat_map, A = accumarray([(WE+1) (SE+1)], 1, [101 101]); end
    else
        if strcmp(which_bout_definition,'standard')
            tempSb = sleep_bouts; tempWb = wake_bouts;
            display(['Number of sleep bouts: ' num2str(length(sleep_bouts))]);
            display(['Mean sleep bout: ' num2str(mean(sleep_bouts))]);
            display(['Number of wake bouts: ' num2str(length(wake_bouts))]);
            display(['Mean sleep bout: ' num2str(mean(wake_bouts))]);
        elseif strcmp(which_bout_definition,'ignore noisy bouts')
            bout_definition = [which_bout_definition, '<' num2str(max_noisy_bout_length) 'sec'];
            [tempSb,tempWb] = calculate_bouts_rnb(nr_events,t,WE,SE,max_noisy_bout_length);
        elseif strcmp(which_bout_definition,'transition path based')
            bout_definition = which_bout_definition;
            [tempSb,tempWb,tempTb] = interval_to_bout(duration_matrix,which_tpb_method);
        end
    end
    
    min_bout = 10;
    nr_of_means = 1;
    
    bs1 = tempSb(tempSb<nr_of_means*mean(tempSb) & tempSb>=min_bout);
    bw1= tempWb(tempWb<nr_of_means*mean(tempWb) & tempWb>=min_bout);
    
    bs2 = tempSb(tempSb>nr_of_means*mean(tempSb) & tempSb<max_bout);
    bw2= tempWb(tempWb>nr_of_means*mean(tempWb) & tempWb<max_bout);
    
    [bsn1,bsx1] = hist(bs1,round((max(bs1)-min(bs1))/quotient1));
    [bwn1,bwx1] = hist(bw1,round((max(bw1)-min(bw1))/quotient1));
    
    [bsn2,bsx2] = hist(bs2,round((max(bs2)-min(bs2))/quotient2));
    [bwn2,bwx2] = hist(bw2,round((max(bw2)-min(bw2))/quotient2));
    
    % --------------------------------------------------------------------
    % -- Estimate sleep and wake scaling parameters
    % Clauset method
    wake_bout_exponent_clauset = (-1)*(1+length(bw1)*1/sum(log(bw1/min_bout)));
    error_wake_exponent = (wake_bout_exponent_clauset-1)/sqrt(length(bw1));
    sleep_bout_exponent_clauset = (-1)*(1+length(bs1)*1/sum(log(bs1/min_bout)));
    error_sleep_exponent = (sleep_bout_exponent_clauset-1)/sqrt(length(bs1));
    display('Clause method');
    display(['Wake bout exponent (of ' num2str(length(bw1)) ' bouts in [' ...
        num2str(min_bout) ',' num2str(nr_of_means*mean(tempWb)) ']) is ' ...
        num2str(wake_bout_exponent_clauset) ' plus/minus ' num2str(error_wake_exponent)]);
    
    display(['Sleep bout exponent (of ' num2str(length(bs1)) ' bouts in [' ...
        num2str(min_bout) ',' num2str(nr_of_means*mean(tempSb)) ']) is ' ...
        num2str(sleep_bout_exponent_clauset) ' plus/minus ' num2str(error_sleep_exponent)]);
    % Linear regression
    md1 = LinearModel.fit(log(bsx1),log(bsn1),'Linear'); coeffs = md1.coefCI;
    sleep_bout_exponent_lin_reg = mean(coeffs(2,:));
    md1 = LinearModel.fit(log(bwx1),log(bwn1),'Linear'); coeffs = md1.coefCI;
    wake_bout_exponent_lin_reg = mean(coeffs(2,:));
    display('Linear regression');
    display(['Wake bout exponent (of ' ...
        num2str(round((max(bw1)-min(bw1))/quotient1)) ' bins in [' ...
        num2str(min_bout) ',' num2str(nr_of_means*mean(tempWb)) ']) is ' ...
        num2str(wake_bout_exponent_lin_reg)]);
    display(['Sleep bout exponent (of ' ...
        num2str(round((max(bs1)-min(bs1))/quotient1)) ' bins in [' ...
        num2str(min_bout) ',' num2str(nr_of_means*mean(tempSb)) ']) is ' ...
        num2str(sleep_bout_exponent_lin_reg)]);
    
    if use_clauset
        exponent_finding_method = 'Clauset algorithm';
        wake_bout_exponent = wake_bout_exponent_clauset;
        sleep_bout_exponent = sleep_bout_exponent_clauset;
    elseif use_linear_regression
        exponent_finding_method = 'Linear regression';
        wake_bout_exponent = wake_bout_exponent_lin_reg;
        sleep_bout_exponent = sleep_bout_exponent_lin_reg;
    end
    % --------------------------------------------------------------------
    
    % Plot bout distributions
    figure,
    subplot(2,2,1); loglog(bsx1,bsn1);
    title(['Short sleep bouts by (' bout_definition '): power-law']);
    x0 = bsx1(1); y0 = bsn1(1); x1 = nr_of_means*mean(tempSb); slope = sleep_bout_exponent;
    hold; plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'k-.'); hold off;
    xlim([0 max_bout]); xlabel(['power law exponent = ' num2str(sleep_bout_exponent) ]);
    
    subplot(2,2,2); semilogy(bsx2,bsn2);
    title(['Long sleep bouts by (' bout_definition '): exponential']);
    x0 = 3e2; y0 = 1e5; x1 = 1.5*1e3; slope = -0.009;
    hold; plot([x0 x1],[y0,y0*exp(slope*(x1-x0))],'k-.'); hold off;
    xlim([0 max_bout]); xlabel(['exponential rate of tail = ' num2str(slope) ]);
    
    subplot(2,2,3); loglog(bwx1,bwn1);
    title(['Short wake bouts by (' bout_definition '): power-law']);
    x0 = bwx1(1); y0 = bwn1(1); x1 = nr_of_means*mean(tempWb); slope = wake_bout_exponent;
    hold; plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'k-.'); hold off;
    xlim([0 max_bout]); xlabel(['power law exponent = ' num2str(wake_bout_exponent) ]);
    
    subplot(2,2,4); semilogy(bwx2,bwn2);
    title(['Long wake bouts by (' bout_definition '): exponential']);
    x0 = 3e2; y0 = 1e5; x1 = 2.5*1e3; slope = -0.004;
    hold; plot([x0 x1],[y0,y0*exp(slope*(x1-x0))],'k-.'); hold off;
    xlim([0 max_bout]); xlabel(['exponential rate of tail = ' num2str(slope) ]);
    
    legend(experiments(which_experiment).info,'Location','Best');
    
    % Plot heat map
    if plot_heat_map
        N = experiments(which_experiment).N;
        wad = experiments(which_experiment).wake_active_domain + 1;
        sad = experiments(which_experiment).sleep_active_domain + 1;
        
        [wpWE,wpSE,~] = find(A==max(max(A(:,1:N/5))));
        [spWE,spSE,~] = find(A==max(max(A(1:N/5,:))));
        
        % 3D
        figure, surf(A); hold;
        plot3(wad(:,2),wad(:,1),A(sub2ind([N+1 N+1],wad(:,1),wad(:,2))),'k*',wpSE,wpWE,A(sub2ind([N+1 N+1],wpWE,wpSE)),'b*');
        plot3(sad(:,2),sad(:,1),A(sub2ind([N+1 N+1],sad(:,1),sad(:,2))),'k*',spSE,spWE,A(sub2ind([N+1 N+1],spWE,spSE)),'b*');
        xlabel('SE'); ylabel('WE'); title(sim_title);
        % 2D
        figure, contour(A,50); hold,
        plot(wad(:,2),wad(:,1),'k*');
        plot(sad(:,2),sad(:,1),'k*');
        xlabel('SE'); ylabel('WE');title(sim_title);
    end
    toc
    
    if plot_full
        quotient = 10;
        [bsn,bsx] = hist(tempSb,round((max(tempSb)-min(tempSb))/quotient));
        [bwn,bwx] = hist(tempWb,round((max(tempWb)-min(tempWb))/quotient));
        
        % Sleep bouts
        figure
        subplot(1,2,1), loglog(bsx,bsn); hold on;
        title(sim_title);
        x0 = min_bout; y0 = bsn(1); x1 = nr_of_means*mean(tempSb);
        slope = sleep_bout_exponent_clauset;
        plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'k-.');% hold off;
        x0 = min_bout; y0 = bsn(1); x1 = nr_of_means*mean(tempSb);
        slope = sleep_bout_exponent_lin_reg;
        plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'k'); hold off;
        legend('Sleep bout pdf',...
            ['exponent by Clauset = ' num2str(sleep_bout_exponent_clauset)],...
            ['exponent by lin. reg. = ' num2str(sleep_bout_exponent_lin_reg)],...
            'Location','Best');
        
        % Wake bouts
        subplot(1,2,2), loglog(bwx,bwn); hold on;
        x0 = min_bout; y0 = bwn(1); x1 = nr_of_means*mean(tempWb);
        slope = wake_bout_exponent_clauset;
        plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'k-.');
        x0 = min_bout; y0 = bwn(1); x1 = nr_of_means*mean(tempWb);
        slope = wake_bout_exponent_lin_reg;
        plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'k');
        legend('Wake bout pdf',...
            ['exponent by Clauset = ' num2str(wake_bout_exponent_clauset)],...
            ['exponent by lin. reg. = ' num2str(wake_bout_exponent_lin_reg)],...
            'Location','Best');
    end
    
elseif compare_all_of_current_round_experiments
    %     % compare_all second round experiments (bimodal heat maps!!!)
    
    plot_heat_maps = 0;
    quotient1 = 5;
    quotient2 = 200;
    min_bout = 0.01;
    nr_of_means = 1;
    max_bout = inf;
    
    for i=1:length(experiments)
        if strcmp(which_bout_definition,'standard')
            bout_definition = which_bout_definition;
            tempSb{i} = experiments(i).sleep_bouts; tempWb{i} = experiments(i).wake_bouts;
        elseif strcmp(which_bout_definition,'transition path based')
            bout_definition = 'tpb - trans ints deleted';
            display('-------------------');
            display(experiments(i).info);
            [tempSb{i},tempWb{i},~] = interval_to_bout(experiments(i).duration_matrix,'intervals are bouts',experiments(i).t(end));
        end
        
        % compare wake and sleep bouts between numerical results
        bs1{i} = tempSb{i}(tempSb{i}<nr_of_means*mean(tempSb{i}) & tempSb{i}>min_bout);
        [Ns1{i},Xs1{i}]=hist(bs1{i},round((max(bs1{i})-min(bs1{i}))/quotient1));
        bw1{i} = tempWb{i}(tempWb{i}<nr_of_means*mean(tempWb{i}) & tempWb{i}>min_bout);
        [Nw1{i},Xw1{i}]=hist(bw1{i},round((max(bw1{i})-min(bw1{i}))/quotient1));
        
        bs2{i} = tempSb{i}(tempSb{i}>nr_of_means*mean(tempSb{i}) & tempSb{i}<max_bout);
        [Ns2{i},Xs2{i}]=hist(bs2{i},round((max(bs2{i})-min(bs2{i}))/quotient2));
        bw2{i} = tempWb{i}(tempWb{i}>nr_of_means*mean(tempWb{i}) & tempWb{i}<max_bout);
        [Nw2{i},Xw2{i}]=hist(bw2{i},round((max(bw2{i})-min(bw2{i}))/quotient2));
        
        % Linear regression analysis for exponent
        md1 = LinearModel.fit(log(Xs1{i}),log(Ns1{i}),'Linear'); coeffs = md1.coefCI;
        sleep_bout_exponent{i} = mean(coeffs(2,:));
        md1 = LinearModel.fit(log(Xw1{i}),log(Nw1{i}),'Linear'); coeffs = md1.coefCI;
        wake_bout_exponent{i} = mean(coeffs(2,:));
        display(['Experiment #' num2str(i)]);
        display(['Wake bout exponent (of ' ...
            num2str(round(max(bw1{i})-min(bw1{i})/quotient1)) ' bins in [' ...
            num2str(min_bout) ',' num2str(nr_of_means*mean(tempWb{i})) ']) is ' ...
            num2str(wake_bout_exponent{i})]);
        display(['Sleep bout exponent (of ' ...
            num2str(round(max(bs1{i})-min(bs1{i})/quotient1)) ' nbins in [' ...
            num2str(min_bout) ',' num2str(nr_of_means*mean(tempSb{i})) ']) is ' ...
            num2str(sleep_bout_exponent{i})]);
    end
    
    figure,
    subplot(2,2,2), semilogy(Xw2{1},Nw2{1},Xw2{2},Nw2{2},Xw2{3},Nw2{3},Xw2{4},Nw2{4}); %,Xw{5},Nw{5});
    title(['Hist of WAKE bouts, ' bout_definition ', semi-log']);
    %     slope=-0.5; hold; plot([min_bout 10*min_bout],[1e3,1e3*(10*min_bout)^(slope)],'k-.'); hold off;
    %     xlim([min_bout 1e4]); xlabel(['Slope of black dashed line = ' num2str(slope) ]);
    legend(experiments(1).info,experiments(2).info,experiments(3).info,experiments(4).info,'Location','Best');
    
    subplot(2,2,4), semilogy(Xs2{1},Ns2{1},Xs2{2},Ns2{2},Xs2{3},Ns2{3},Xs2{4},Ns2{4}); %,Xs{5},Ns{5});
    title(['Hist of SLEEP bouts, ' bout_definition ', semi-log']);
    %     slope=-0.5; hold; plot([min_bout 10*min_bout],[1e3,1e3*(10*min_bout)^(slope)],'k-.'); hold off;
    %     xlim([min_bout 1e4]); xlabel(['Slope of black dashed line = ' num2str(slope) ]);
    subplot(2,2,1), loglog(Xw1{1},Nw1{1},Xw1{2},Nw1{2},Xw1{3},Nw1{3},Xw1{4},Nw1{4}); %,Xw{5},Nw{5});
    title(['Hist of WAKE bouts, ' bout_definition ', log-log']);
    x0 = 2; y0 = 5e2; x1 = 5e2; slope = -0.8;
    hold; plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'k-.'); hold off;
    xlim([0 max_bout]); xlabel(['power law exponent = ' num2str(slope) ]);
    subplot(2,2,3), loglog(Xs1{1},Ns1{1},Xs1{2},Ns1{2},Xs1{3},Ns1{3},Xs1{4},Ns1{4}); %,Xs{5},Ns{5});
    title(['Hist of SLEEP bouts, ' bout_definition ', log-log']);
    x0 = 2; y0 = 5e2; x1 = 5e2; slope = -0.8;
    hold; plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'k-.'); hold off;
    xlim([0 max_bout]); xlabel(['power law exponent = ' num2str(slope) ]);
    
    if plot_heat_maps
        figure,
        for k=1:length(experiments)
            subplot(2,2,k); surf(experiments(k).A);
            xlim([0 experiments(k).N+1]); ylim([0 experiments(k).N+1]);
            if k==1, title('Blue'); elseif k==2, title('Green'); elseif k==3, title('Red'); elseif k==4, title('Light Blue'); end
        end
    end
    
    figure,
    loglog(Xw1{1},Nw1{1},Xw1{2},Nw1{2},Xw1{3},Nw1{3},Xw1{4},Nw1{4});
    title(['Hist of WAKE bouts, ' bout_definition ', log-log']);
    legend([num2str(wake_bout_exponent{1}) ' [' experiments(1).info ']'],...
        [num2str(wake_bout_exponent{2}) ' [' experiments(2).info ']'],...
        [num2str(wake_bout_exponent{3}) ' [' experiments(3).info ']'],...
        [num2str(wake_bout_exponent{4}) ' [' experiments(4).info ']'],...
        'Location','Best');
    
    figure,
    loglog(Xs1{1},Ns1{1},Xs1{2},Ns1{2},Xs1{3},Ns1{3},Xs1{4},Ns1{4});
    title(['Hist of SLEEP bouts, ' bout_definition ', log-log']);
    legend([num2str(sleep_bout_exponent{1}) ' [' experiments(1).info ']'],...
        [num2str(sleep_bout_exponent{2}) ' [' experiments(2).info ']'],...
        [num2str(sleep_bout_exponent{3}) ' [' experiments(3).info ']'],...
        [num2str(sleep_bout_exponent{4}) ' [' experiments(4).info ']'],...
        'Location','Best');
    
elseif plot_sample_trajectories_during_bouts
    which_experiment = 2;
    make_movie_from_bout_trajs = 0;
    what_type_of_bouts = [0];
    
    min_bout_size = 0;
    max_bout_size = 1000;
    nr_bout_trajs = 50;
    pause_after_figure = 1;
    
    if make_movie_from_bout_trajs, close all, end
    
    wake_domain = experiments(which_experiment).wake_active_domain;
    sleep_domain = experiments(which_experiment).sleep_active_domain;
    
    for i=1:nr_bout_trajs
        bout_index_set = find(experiments(which_experiment).duration_matrix(:,1) > ...
            min_bout_size & experiments(which_experiment).duration_matrix(:,1) < ...
            max_bout_size);
        bout_index = bout_index_set(randperm(length(bout_index_set),1));
        bout_type = experiments(which_experiment).duration_matrix(bout_index,2);
        while ~ismember(bout_type,what_type_of_bouts)
            bout_index = bout_index_set(randperm(length(bout_index_set),1));
            bout_type = experiments(which_experiment).duration_matrix(bout_index,2);
        end
        bout_begin_time = sum(experiments(which_experiment).duration_matrix(1:bout_index-1,1));
        bout_end_time = bout_begin_time + ...
            experiments(which_experiment).duration_matrix(bout_index,1);
        
        traj_index_set_min = find(experiments(which_experiment).t >= ...
            bout_begin_time,1);
        traj_index_set_max = find(experiments(which_experiment).t > ...
            bout_end_time,1);
        traj_index_set = (traj_index_set_min:traj_index_set_max-1)';
        
        % Plot on 2D: WE vs SE
        
        figure,  hold;
        plot(wake_domain(:,2),wake_domain(:,1),'cx');
        plot(sleep_domain(:,2),sleep_domain(:,1),'cx');
        plot(experiments(which_experiment).SE(traj_index_set),...
            experiments(which_experiment).WE(traj_index_set),'b.-');
        plot(experiments(which_experiment).SE(traj_index_set(1)),...
            experiments(which_experiment).WE(traj_index_set(1)),'r*');
        plot(experiments(which_experiment).SE(traj_index_set(end)),...
            experiments(which_experiment).WE(traj_index_set(end)),'k*');
        axis([0 101 0 101]); xlabel('SE'); ylabel('WE');
        legend(experiments(which_experiment).info,'Location','Best');
        if bout_type==-1
            str_bout_type = 'sleep';
        elseif bout_type==1
            str_bout_type = 'wake';
        elseif bout_type==0
            str_bout_type = 'transition';
        end
        title(['Sample ' str_bout_type ' bout, size = ' ...
            num2str(experiments(which_experiment).duration_matrix(bout_index,1))]);
        
        if pause_after_figure, pause; end
    end
    if make_movie_from_bout_trajs
        vidObj = VideoWriter(['bout_trajs_' datestr(now,'mmmdd_HHMM') '.avi']);
        vidObj.FrameRate = 0.5; % 2 seconds per frame
        open(vidObj);
        for k=1:nr_bout_trajs, writeVideo(vidObj,getframe(figure(k))); end
        close(vidObj);
    end
elseif compare_two_experiments
    % compare_experiments i and j
    if ~use_current_workspace
        expi = experiments(i);
        expj = experiments(j);
        
        cluster1 = expi.info;
        cluster2 = expj.info;
        
        if strcmp(which_bout_definition,'standard')
            bout_definition = which_bout_definition;
            tempSb1 = expi.sleep_bouts; tempWb1 = expi.wake_bouts;
            tempSb2 = expj.sleep_bouts; tempWb2 = expj.wake_bouts;
        elseif strcmp(which_bout_definition,'threshold')
            bout_definition = [which_bout_definition, '=' num2str(threshold)];
            [tempSb1,tempWb1] = calculate_bouts_threshold(expi.nr_events,expi.t,expi.WE,expi.SE,expi.N,threshold);
            [tempSb2,tempWb2] = calculate_bouts_threshold(expj.nr_events,expj.t,expj.WE,expj.SE,expj.N,threshold);
        elseif strcmp(which_bout_definition,'ignore noisy bouts')
            bout_definition = [which_bout_definition, '<' num2str(max_noisy_bout_length) 'sec'];
            [tempSb1,tempWb1] = calculate_bouts_rnb(expi.nr_events,expi.t,expi.WE,expi.SE,max_noisy_bout_length);
            [tempSb2,tempWb2] = calculate_bouts_rnb(expj.nr_events,expj.t,expj.WE,expj.SE,max_noisy_bout_length);
        elseif strcmp(which_bout_definition,'compare two definitions')
            bout_definition = which_bout_definition;
            tempSb1 = expi.sleep_bouts; tempWb1 = expi.wake_bouts;
            [tempSb2,tempWb2] = calculate_bouts_rnb(expj.nr_events,expj.t,expj.WE,expj.SE,max_noisy_bout_length);
        elseif strcmp(which_bout_definition,'transition path based')
            bout_definition = 'tpb - trans ints deleted';
            [tempSb1,tempWb1,~] = interval_to_bout(expi.duration_matrix,'intervals are bouts',expi.t(end));
            [tempSb2,tempWb2,~] = interval_to_bout(expj.duration_matrix,'intervals are bouts',expj.t(end));
        end
    else
        bout_definition = which_bout_definition;
        cluster1 = ['[standard]:' sim_title];
        cluster2 = ['[tpb]:' sim_title];
        tempSb1 = sleep_bouts; tempWb1 = wake_bouts;
        [tempSb2,tempWb2,tempTb2] = interval_to_bout(duration_matrix,which_tpb_method,t(end));
    end
    
    % compare wake and sleep bouts between numerical results
    bs1 = tempSb1(tempSb1<nr_of_means*mean(tempSb1) & tempSb1>m);
    bs2 = tempSb2(tempSb2<nr_of_means*mean(tempSb2) & tempSb2>m);
    [Ns1,Xs1]=hist(bs1,round(max(bs1)-min(bs1))/quotient);
    [Ns2,Xs2]=hist(bs2,round(max(bs2)-min(bs2))/quotient);
    
    bw1 = tempWb1(tempWb1<nr_of_means*mean(tempWb1) & tempWb1>m);
    bw2 = tempWb2(tempWb2<nr_of_means*mean(tempWb2) & tempWb2>m);
    [Nw1,Xw1]=hist(bw1,round(max(bw1)-min(bw1))/quotient);
    [Nw2,Xw2]=hist(bw2,round(max(bw2)-min(bw2))/quotient);
    
    figure,
    subplot(2,2,2), semilogy(Xw1,Nw1,'b-',Xw2,Nw2,'r-');
    title(['Hist of WAKE bouts, ' bout_definition ', semi-log']);
    legend(cluster1,cluster2,'Location','Best');
    subplot(2,2,4), semilogy(Xs1,Ns1,'b-',Xs2,Ns2,'r-');
    title(['Hist of SLEEP bouts, ' bout_definition ', semi-log']);
    
    subplot(2,2,1), loglog(Xw1,Nw1,'b-',Xw2,Nw2,'r-');
    title(['Hist of WAKE bouts, ' bout_definition ', log-log']);
    subplot(2,2,3), loglog(Xs1,Ns1,'b-',Xs2,Ns2,'r-');
    title(['Hist of SLEEP bouts, ' bout_definition ', log-log']);
    
    % Plot sample WE vs SE states for the two experiments in the same figure window
    if plot_sample_simulations
        figure,
        count = 1;
        for k=[i,j]
            
            a = experiments(k);
            t = a.t; WE = a.WE; SE = a.SE; WI = a.WI; SI = a.SI; N = a.N;
            t_out = a.t_out; y_out = a.y_out; dW = a.dW; dI = a.dI; dS = a.dS;
            wake_cluster = a.wake_cluster; sleep_cluster = a.sleep_cluster;
            
            t_MF = find(t_out<max(t(X)));
            subplot(2,1,count), plot(t(X),WE(X)/N,'b',t_out(t_MF),y_out(t_MF,2),'r',...
                t(X),SE(X)/N,'g',t_out(t_MF),y_out(t_MF,4),'r'); %,...
            title(['W=' wake_cluster '(' num2str(N) ',' num2str(dW) ')'...
                ', S=' sleep_cluster '(' num2str(N) ',' num2str(dS) ')'...
                ', dI=' num2str(dI)]);
            legend('Wake-Excited','MF-eq','Sleep-Excited','MF-eq','Location','Best');
            count = count+1;
            
        end
    end
elseif compare_tpb_methods_and_standard
    % Calculate bouts from intervals using three different ways
    [tempSb1,tempWb1,tempTb1] = interval_to_bout(duration_matrix,'intervals are bouts',t(end));
    [tempSb2,tempWb2,tempTb2] = interval_to_bout(duration_matrix,'transition intervals are added to previous bout',t(end));
    [tempSb3,tempWb3,tempTb3] = interval_to_bout(duration_matrix,'transition intervals are cut in half and added to closest wake/sleep interval',t(end));
    
    % Bouts calculated by the standard definition
    tempSb4 = sleep_bouts; tempWb4 = wake_bouts;
    display(['Number of sleep bouts: ' num2str(length(tempSb4))]);
    display(['Number of wake bouts: ' num2str(length(tempWb4))]);
    display(['Number of transition bouts: ' num2str(length([]))]);
    
    display(['Mean sleep bout: ' num2str(mean(tempSb4))]);
    display(['Mean wake bout: ' num2str(mean(tempWb4))]);
    display(['Mean transition bout: ' num2str(mean([]))]);
    
    display(['Fraction of sleep time: ' num2str(sum(tempSb4)/t(end))]);
    display(['Fraction of wake time: ' num2str(sum(tempWb4)/t(end))]);
    display(['Fraction of transition time: ' num2str(sum([])/t(end))]);
    
    % Compare the bout distributions for the above three different ways
    bs1 = tempSb1(tempSb1<nr_of_means*mean(tempSb1) & tempSb1>m);
    bs2 = tempSb2(tempSb2<nr_of_means*mean(tempSb2) & tempSb2>m);
    bs3 = tempSb3(tempSb3<nr_of_means*mean(tempSb3) & tempSb3>m);
    bs4 = tempSb4(tempSb4<nr_of_means*mean(tempSb4) & tempSb4>m);
    
    [Ns1,Xs1]=hist(bs1,round(max(bs1)-min(bs1))/quotient);
    [Ns2,Xs2]=hist(bs2,round(max(bs2)-min(bs2))/quotient);
    [Ns3,Xs3]=hist(bs3,round(max(bs3)-min(bs3))/quotient);
    [Ns4,Xs4]=hist(bs4,round(max(bs4)-min(bs4))/quotient);
    
    bw1 = tempWb1(tempWb1<nr_of_means*mean(tempWb1) & tempWb1>m);
    bw2 = tempWb2(tempWb2<nr_of_means*mean(tempWb2) & tempWb2>m);
    bw3 = tempWb3(tempWb3<nr_of_means*mean(tempWb3) & tempWb3>m);
    bw4 = tempWb4(tempWb4<nr_of_means*mean(tempWb4) & tempWb4>m);
    
    [Nw1,Xw1]=hist(bw1,round(max(bw1)-min(bw1))/quotient);
    [Nw2,Xw2]=hist(bw2,round(max(bw2)-min(bw2))/quotient);
    [Nw3,Xw3]=hist(bw3,round(max(bw3)-min(bw3))/quotient);
    [Nw4,Xw4]=hist(bw4,round(max(bw4)-min(bw4))/quotient);
    
    figure,
    subplot(2,2,2), semilogy(Xw1,Nw1,'b-',Xw2,Nw2,'r-',Xw3,Nw3,'-g',Xw4,Nw4,'-k');
    title('Hist of WAKE bouts, semi-log');
    legend('ints=bouts','ints added to prev','ints are cut in half','standard','Location','Best');
    subplot(2,2,4), semilogy(Xs1,Ns1,'b-',Xs2,Ns2,'r-',Xs3,Ns3,'-g',Xs4,Ns4,'-k');
    title('Hist of SLEEP bouts, semi-log');
    subplot(2,2,1), loglog(Xw1,Nw1,'b-',Xw2,Nw2,'r-',Xw3,Nw3,'-g',Xw4,Nw4,'-k');
    title('Hist of WAKE bouts, log-log');
    subplot(2,2,3), loglog(Xs1,Ns1,'b-',Xs2,Ns2,'r-',Xs3,Ns3,'-g',Xs4,Ns4,'-k');
    title('Hist of SLEEP bouts, log-log');
    
elseif plot_drift_diffusion_peclet
    %%
    which_experiment = 1;
    
    k = 100; % k^2 is the number of points reprensented in drift and diffusion plots
    
    plot_sample_trajectories = 0;
    plot_null_clines = 0;
    
    if ~use_current_workspace
        N = experiments(which_experiment).N;
        sim_title = experiments(which_experiment).info;
        F = experiments(which_experiment).F;
        D = experiments(which_experiment).D;
    end
    % Plot State vs Truncated norm of drift
    normF = zeros(N+1);
    for i=1:N+1
        for j=1:N+1
            normF(i,j) = norm([F(i,j,1) F(i,j,2)]);
        end
    end
    normF_truncated = normF;
    normF_truncated(normF_truncated>1)=1;
    figure,  hold;
    [c, h] = contour(1:N+1,1:N+1,normF_truncated);
    clabel(c, h);
    title(['State vs 2-norm of Drift for' sim_title]);
    xlabel('SE'); ylabel('WE');
    
    % Add sample trajectories
    if plot_sample_trajectories
        for i=0.1:0.1:0.9
            for j=0.1:0.1:0.9
                [~,y_out] = MF_model(0.1,i,0.1,j,...
                    experiments(which_experiment).dW,...
                    experiments(which_experiment).dS,...
                    experiments(which_experiment).dIW,...
                    experiments(which_experiment).dIS,...
                    0.001,0.003,0.016,0.002,0.005);
                plot(100*y_out(:,4),100*y_out(:,2)); % plot SE vs WE
            end
        end
    end
    
    % Add MF-equation null clines data from xppaut
    if plot_null_clines
        fid = fopen('WE_nc.txt','rt');
        A = textscan(fid,'%f %f %d','HeaderLines',1);
        WEncx = A{1}; WEncy = A{2};
        plot(WEncx*100,WEncy*100,'r','Linewidth',2);
        fclose(fid);
        fid = fopen('SE_nc.txt','rt');
        A = textscan(fid,'%f %f %d','HeaderLines',1);
        SEncx = A{1}; SEncy = A{2};
        plot(SEncx*100,SEncy*100,'r','Linewidth',2);
        fclose(fid);
    end
    
    U = F(:,:,1,1); V = F(:,:,2,1);
    U_new=zeros(k,k); V_new=zeros(k,k);
    for i=1:k
        for j=1:k
            U_new(i,j)=sum(sum(U(1+(i-1)*N/k:i*N/k,1+(j-1)*N/k:j*N/k)));
            V_new(i,j)=sum(sum(V(1+(i-1)*N/k:i*N/k,1+(j-1)*N/k:j*N/k)));
        end
    end
    [X,Y] = meshgrid(0.5*N/k:N/k:N,0.5*N/k:N/k:N);
    figure, quiver(X,Y,U_new,V_new); axis([0 N 0 N]);
    title(['Drift vector field (' sim_title ')']);
    xlabel('SE'); ylabel('WE');
    
    % Plotting diffusion tensor
    eigV = zeros(2,2,k,k); % 4D - eigenvector array
    eigD = zeros(2,2,k,k); % 4D - eigenvalue array
    
    scaling=1;
    figure, axis([0 N 0 N]);
    title(['Diffusion tensors for (' sim_title ')']);
    xlabel('SE'); ylabel('WE');
    hold on;
    
    for i=1:k
        for j=1:k
            [eigV(:,:,i,j),eigD(:,:,i,j)] = ...
                eig(mean(mean(D(:,:,1+(i-1)*N/k:i*N/k,1+(j-1)*N/k:j*N/k),3),4));
            
            vv = eigV(:,:,i,j); % eigenvectors: vv(:,1) and vv(:,2)
            dd = eigD(:,:,i,j); % eigenvalues: dd(1,1) and dd(2,2)
            plot([(i-0.5)*N/k-scaling*dd(1,1)/2*vv(1,1),...
                (i-0.5)*N/k+scaling*dd(1,1)/2*vv(1,1)],...
                [(j-0.5)*N/k-scaling*dd(1,1)/2*vv(2,1),...
                (j-0.5)*N/k+scaling*dd(1,1)/2*vv(2,1)]);
            plot([(i-0.5)*N/k-scaling*dd(2,2)/2*vv(1,2),...
                (i-0.5)*N/k+scaling*dd(2,2)/2*vv(1,2)],...
                [(j-0.5)*N/k-scaling*dd(2,2)/2*vv(2,2),...
                (j-0.5)*N/k+scaling*dd(2,2)/2*vv(2,2)]);
            
            
        end
    end
    hold off;
    
    % Plot Peclet number for 2D system
    peclet = zeros(k);
    for i=1:k
        for j=1:k
            dd = eigD(:,:,j,i); % eigenvalues: dd(1,1) and dd(2,2)
            peclet(i,j) = N*norm([U_new(i,j) V_new(i,j)])/max(dd(1,1),dd(2,2));
        end
    end
    peclet_truncated = peclet; peclet_truncated(peclet>50) = 50;
    [c, h] = contour(peclet_truncated);
    clabel(c, h);
    title(['Peclet number for 2D drift-diffusion estimation of' sim_title]);
    xlabel('SE'); ylabel('WE');
    
    
elseif compare_transition_intervals
    quotient1 = 1;
    quotient2 = 4;
    nr_of_means = 1;
    
    for i = [1 2 4 5] %1:length(experiments)
        if strcmp(which_bout_definition,'standard')
            error('Trying to plot transition intervals and bout defn is standard!!!');
        elseif strcmp(which_bout_definition,'transition path based')
            bout_definition = 'tpb - trans ints deleted';
            [~,~,trans_int{i}] = interval_to_bout(...
                experiments(i).duration_matrix,'intervals are bouts',...
                experiments(i).t(end));
        end
        [nt1{i},xt1{i}]=hist(trans_int{i},...
            round(max(trans_int{i})-min(trans_int{i}))/quotient1);
        ti2{i} = trans_int{i}(trans_int{i}>nr_of_means*mean(trans_int{i}));
        [nt2{i},xt2{i}]=hist(ti2{i},round(max(ti2{i})-min(ti2{i}))/quotient2);
    end
    figure,
    subplot(2,1,1); plot(xt1{1},nt1{1},xt1{2},nt1{2},xt1{4},nt1{4},xt1{5},nt1{5});
    title('Distribution of transition intervals');
    legend(experiments(1).info,experiments(2).info,...
        experiments(4).info,experiments(5).info,'Location','Best');
    subplot(2,1,2); semilogy(xt2{1},nt2{1},xt2{2},nt2{2},xt2{4},nt2{4},xt2{5},nt2{5});
    title('Exponential tail in transition intervals');
    hold off;
elseif compare_statistics
    plot_ws_sw_transitions = 0;
    
    for i=1:length(experiments)
        [sb{i},wb{i},tb{i},wstb{i},swtb{i}] = interval_to_bout(...
            experiments(i).duration_matrix,which_tpb_method,...
            experiments(i).t(end));
        mean_sw_cycle(i) = experiments(i).t(end)/length(sb{i});
        sb_frac(i) = sum(sb{i})/experiments(i).t(end);
        wb_frac(i) = sum(wb{i})/experiments(i).t(end);
        tb_frac(i) = sum(tb{i})/experiments(i).t(end);
        
        tb_mean(i) = mean(tb{i});
        
        if plot_ws_sw_transitions
            [nwstb,xwstb] = hist(wstb{i},20);
            [nswtb,xswtb] = hist(swtb{i},20);
            
            figure, plot(xwstb,nwstb,'b*-',xswtb,nswtb,'r*-');
            legend('W->S transition histogram','S->W transition histogram',...
                'Location','Best');
            title(experiments(i).info);
        end
    end
    
    sampletitle = ['ER-ER(' num2str(experiments(1).N) ',' ...
        num2str(experiments(1).dS) ',' num2str(experiments(1).dIS) ')'];
    
    figure, subplot(2,2,1)
    axis([0 100 0 100]); hold on;
    plot(sb_frac(1:3:end)*100,wb_frac(1:3:end)*100,'b*');
    plot(sb_frac(2:3:end)*100,wb_frac(2:3:end)*100,'g*');
    plot(sb_frac(3:3:end)*100,wb_frac(3:3:end)*100,'r*');
    xlabel('Fraction of total sleep time');
    ylabel('Fraction of total wake time'); hold off;
    
    subplot(2,2,2), hold on;
    plot(tb_mean(1:3:end),mean_sw_cycle(1:3:end),'b*');
    plot(tb_mean(2:3:end),mean_sw_cycle(2:3:end),'g*');
    plot(tb_mean(3:3:end),mean_sw_cycle(3:3:end),'r*');
    legend(sampletitle,'(wake) ER-SF (sleep)','SF-SF','Location','Best');
    xlabel('Mean transition time');
    ylabel('Mean sleep-wake cycle time'); hold off;
    
    
    subplot(2,2,3), xlim([0 100]); hold on;
    plot(sb_frac(1:3:end)*100,mean_sw_cycle(1:3:end),'b*');
    plot(sb_frac(2:3:end)*100,mean_sw_cycle(2:3:end),'g*');
    plot(sb_frac(3:3:end)*100,mean_sw_cycle(3:3:end),'r*');
    xlabel('Fraction of total sleep time');
    ylabel('Mean sleep-wake cycle time'); hold off;
    
    subplot(2,2,4), xlim([0 100]); hold on;
    plot(wb_frac(1:3:end)*100,mean_sw_cycle(1:3:end),'b*');
    plot(wb_frac(2:3:end)*100,mean_sw_cycle(2:3:end),'g*');
    plot(wb_frac(3:3:end)*100,mean_sw_cycle(3:3:end),'r*');
    xlabel('Fraction of total wake time');
    ylabel('Mean sleep-wake cycle time'); hold off;
    
    
elseif compute_SW_WS_transition_probability_densities
    % Start with a network configuration and predetermined wake and sleep
    % active domains. Initialize the process and run it for a while in
    % order to equilibrate. Next run the process from this equilibrium for
    % many times to get a distribution of P(activation switch in less than
    % T time units). Do this for ER and SF and compare the two densities
    
    % which_experiment = 7;
    % Use ER-SF to see difference in SF and ER, hence use 2+3k
    save_variables = 1;
    
    a = experiments(which_experiment);
    display(a.info);
    
    end_events = (1e2:1e2:3e3)';
    nr_trials = 10;
    nr_initial_state_realizations = 10;
    WS_success_rate_matrix = zeros(length(end_events),2,nr_initial_state_realizations);
    SW_success_rate_matrix = zeros(length(end_events),2,nr_initial_state_realizations);
    
    tic
    for k = 1:nr_initial_state_realizations
        % Approximate pdf of W->S transition times
        WS_success_rate_matrix(:,:,k) = two_cluster_simulation_until_switch(a.dW,a.WW,a.dS,a.SS,...
            a.dIW,a.WS,a.dIS,a.SW,a.wake_active_domain,a.sleep_active_domain,...
            'W->S',end_events,nr_trials);
        
        % Approximate pdf of S->W transition times
        SW_success_rate_matrix(:,:,k) = two_cluster_simulation_until_switch(a.dW,a.WW,a.dS,a.SS,...
            a.dIW,a.WS,a.dIS,a.SW,a.wake_active_domain,a.sleep_active_domain,...
            'S->W',end_events,nr_trials);
    end
    toc
    
    figure, hold on;
    for i=1:nr_initial_state_realizations
        plot(WS_success_rate_matrix(:,1,i),WS_success_rate_matrix(:,2,i),'b.-');
        plot(SW_success_rate_matrix(:,1,i),SW_success_rate_matrix(:,2,i),'r.-');
    end
    legend('W->S transition cdf','S->W transition cdf','Location','Best');
    xlabel(['Sim# ' num2str(which_experiment) ': ' a.info]);
    hold off;
    
    % Plot W->S and S->W transition cdf's averaged over different initial
    % states in home domains
    mean_WS_srm = mean(WS_success_rate_matrix,3);
    mean_SW_srm = mean(SW_success_rate_matrix,3);
    figure, hold on;
    plot(mean_WS_srm(:,1),mean_WS_srm(:,2),'b*-');
    plot(mean_SW_srm(:,1),mean_SW_srm(:,2),'r*-');
    legend('W->S transition cdf averaged over initial states',...
        'S->W transition cdf averaged over initial states','Location','Best');
    title(['CDFs Averaged over ' num2str(nr_initial_state_realizations) ...
        ' initial state']);
    xlabel(['Sim# ' num2str(which_experiment) ': ' a.info]);
    hold off;
    
    
    [sb,wb,tb,wstb,swtb] = interval_to_bout(a.duration_matrix,...
        which_tpb_method,a.t(end));
    
    % Plot P(S->W transition|S) divided by P(W->S transition|W)
    figure, hold on;
    plot(end_events,mean_WS_srm(:,2)./mean_SW_srm(:,2),'b*-');
    plot(end_events,ones(length(end_events),1)*mean(sb)/mean(wb),'k.-.');
    legend('P(W->S transition|W) divided by P(S->W transition in |S)',...
        'Mean sleep bout divided by mean wake bout','Location','Best');
    xlabel(['Sim# ' num2str(which_experiment) ': ' a.info]);
    hold off;
    
    if save_variables
        filename = ['SW_WS_trans_cdfs_' datestr(now,'mmmdd_HHMM') '.mat'];
        save(filename,'WS_success_rate_matrix','SW_success_rate_matrix',...
            'which_experiment','nr_trials','nr_initial_state_realizations');
    end
    
    
elseif plot_WE_SE_trajectories_vs_t
    which_experiment = 1;
    how_many = 10;
    middle_lines_for_transition_regions = 1;
    
    a = experiments(which_experiment);
    
    new_t = cumsum(a.duration_matrix(:,1));
    for i=0:how_many
        figure
        X1 = i*2e4+1:(i+1)*2e4;
        plot(a.t(X1),a.WE(X1),'r',a.t(X1),a.SE(X1),'b');
        title(a.info); xlabel('t');
        legend('WE','SE','Location','Best');
        

    end
    
elseif plot_sleep_wake_survivals_and_WS_SW_transition_cdfs
    
    which_experiment = 8;
    
    if which_experiment == 5, load SW_WS_trans_cdfs_Sep15_0336.mat,
    elseif which_experiment == 8, load SW_WS_trans_cdfs_Sep15_0559.mat; end
    
    % Plot W->S and S->W transition cdf's averaged over different initial
    % states in home domains
    time_per_event = experiments(which_experiment).t(end)/...
        experiments(which_experiment).nr_events;
    mean_WS_srm = mean(WS_success_rate_matrix,3);
    mean_SW_srm = mean(SW_success_rate_matrix,3);
    figure, hold on;
    plot(mean_WS_srm(:,1)*time_per_event,mean_WS_srm(:,2),'b*-');
    plot(mean_SW_srm(:,1)*time_per_event,mean_SW_srm(:,2),'r*-');
    legend('W->S transition cdf averaged over initial states',...
        'S->W transition cdf averaged over initial states','Location','Best');
    title(['CDFs Averaged over ' num2str(nr_initial_state_realizations) ...
        ' initial state']);
    % xlabel('Number of events');
    xlabel('t: Simulation time'); ylabel('Prob(transition occurs before t)')
    hold off;
    
    % Plot cdf's for wake bout size and sleep bout size
    quotient = 1;
    [sb,wb,~,~,~] = interval_to_bout(experiments(i).duration_matrix,...
        which_tpb_method,experiments(i).t(end));
    
    [sbn,sbx] = hist(sb,round((max(sb)-min(sb))/quotient));
    [wbn,wbx] = hist(wb,round((max(wb)-min(wb))/quotient));
    cum_sbn = cumsum(sbn)/sum(sbn); cum_wbn = cumsum(wbn)/sum(wbn);
    figure, hold on;
    plot(wbx,cum_wbn,'b'); plot(sbx,cum_sbn,'r');
    hold off;
    
elseif adaptive_power_law_detection
    
    which_experiment = 3;
    
    quotient1 = 1;
    quotient2 = 100;
    
    plot_full = 1;
    bout_definition = which_bout_definition;
    which_tpb_method = 'intervals are bouts';%'transition intervals are added to previous bout';
    tic
    if ~use_current_workspace
        sleep_bouts = experiments(which_experiment).sleep_bouts;
        wake_bouts = experiments(which_experiment).wake_bouts;
        
        sim_title = experiments(which_experiment).info;
        A = experiments(which_experiment).A;
        WE = experiments(which_experiment).WE;
        SE = experiments(which_experiment).SE;
        nr_events = experiments(which_experiment).nr_events;
        t = experiments(which_experiment).t;
        
        [tempSb,tempWb,tempTb] = interval_to_bout(...
            experiments(which_experiment).duration_matrix,which_tpb_method);
        
        % if plot_heat_map, A = accumarray([(WE+1) (SE+1)], 1, [101 101]); end
    else
        if strcmp(which_bout_definition,'standard')
            tempSb = sleep_bouts; tempWb = wake_bouts;
            display(['Number of sleep bouts: ' num2str(length(sleep_bouts))]);
            display(['Mean sleep bout: ' num2str(mean(sleep_bouts))]);
            display(['Number of wake bouts: ' num2str(length(wake_bouts))]);
            display(['Mean sleep bout: ' num2str(mean(wake_bouts))]);
        elseif strcmp(which_bout_definition,'ignore noisy bouts')
            bout_definition = [which_bout_definition, '<' num2str(max_noisy_bout_length) 'sec'];
            [tempSb,tempWb] = calculate_bouts_rnb(nr_events,t,WE,SE,max_noisy_bout_length);
        elseif strcmp(which_bout_definition,'transition path based')
            bout_definition = which_bout_definition;
            [tempSb,tempWb,tempTb] = interval_to_bout(duration_matrix,which_tpb_method);
        end
    end
    
    min_nr_of_means = 0.8;
    max_nr_of_means = 1;
    nr_max_bouts = 50;
    
    max_bouts_vector = linspace(min_nr_of_means*mean([mean(tempSb),mean(tempWb)]),...
        max_nr_of_means*mean([mean(tempSb),mean(tempWb)]),...
        nr_max_bouts);
    min_bouts_vector = [50 60 70 80];
    
    wake_bout_exponent_lin_reg = zeros(nr_max_bouts,length(min_bouts_vector));
    sleep_bout_exponent_lin_reg = zeros(nr_max_bouts,length(min_bouts_vector));
    wake_bout_exponent_clauset = zeros(nr_max_bouts,length(min_bouts_vector));
    sleep_bout_exponent_clauset = zeros(nr_max_bouts,length(min_bouts_vector));
    wake_bout_exponent_MLE = zeros(nr_max_bouts,length(min_bouts_vector));
    sleep_bout_exponent_MLE = zeros(nr_max_bouts,length(min_bouts_vector));
    
    for k=1:length(min_bouts_vector)
        min_bout = min_bouts_vector(k);
        for i=1:nr_max_bouts
            max_bout = max_bouts_vector(i);
            bs1 = tempSb(tempSb<max_bout & tempSb>=min_bout);
            bw1= tempWb(tempWb<max_bout & tempWb>=min_bout);
            
            [bsn1,bsx1] = hist(bs1,round((max(bs1)-min(bs1))/quotient1));
            [bwn1,bwx1] = hist(bw1,round((max(bw1)-min(bw1))/quotient1));
            
            % Clauset method
            wake_bout_exponent_clauset(i,k) = (-1)*(1+length(bw1)*1/sum(log(bw1/min_bout)));
            sleep_bout_exponent_clauset(i,k) = (-1)*(1+length(bs1)*1/sum(log(bs1/min_bout)));
            
            % Linear regression
            md1 = LinearModel.fit(log(bsx1),log(bsn1),'Linear'); coeffs = md1.coefCI;
            sleep_bout_exponent_lin_reg(i,k) = mean(coeffs(2,:));
            md1 = LinearModel.fit(log(bwx1),log(bwn1),'Linear'); coeffs = md1.coefCI;
            wake_bout_exponent_lin_reg(i,k) = mean(coeffs(2,:));
            
            % MLE
            % Sleep bout exponent
            set_of_alphas = (0.1:0.01:1.5)';
            U = ((max_bout).^(1-set_of_alphas)-...
                min_bout.^(1-set_of_alphas))./(1-set_of_alphas);
            U(isnan(U)|isinf(U)) = log(max_bout/min_bout);
            L_of_alphas = (-1)*set_of_alphas*sum(log(bs1)) - length(bs1)*log(U);
            [~,L_max_index] = max(L_of_alphas);
            sleep_bout_exponent_MLE(i,k) = (-1)*set_of_alphas(L_max_index);
            % Wake bout exponent
            U = ((max_bout).^(1-set_of_alphas)-...
                min_bout.^(1-set_of_alphas))./(1-set_of_alphas);
            U(isnan(U)|isinf(U)) = log(max_bout/min_bout);
            L_of_alphas = (-1)*set_of_alphas*sum(log(bw1)) - length(bw1)*log(U);
            [~,L_max_index] = max(L_of_alphas);
            wake_bout_exponent_MLE(i,k) = (-1)*set_of_alphas(L_max_index);
        end
    end
    % Plot exponents found by three methods
    figure, hold on;
    for k =1:length(min_bouts_vector)
        plot(max_bouts_vector,sleep_bout_exponent_lin_reg(:,k),'b*-',...
            max_bouts_vector,sleep_bout_exponent_clauset(:,k),'r*-',...
            max_bouts_vector,sleep_bout_exponent_MLE(:,k),'k*-');
        text(max_bouts_vector(end)+1,sleep_bout_exponent_lin_reg(end,k),num2str(min_bouts_vector(k)));
        text(max_bouts_vector(end)+1,sleep_bout_exponent_clauset(end,k),num2str(min_bouts_vector(k)));
        text(max_bouts_vector(end)+1,sleep_bout_exponent_MLE(end,k),num2str(min_bouts_vector(k)));
    end
    title(['Sleep bout exponents for ' sim_title]);
    legend('Linear reg','Clauset','MLE','Location','Best');
    hold off;
    
    figure, hold on;
    for k =1:length(min_bouts_vector)
        plot(max_bouts_vector,wake_bout_exponent_lin_reg(:,k),'b*-',...
            max_bouts_vector,wake_bout_exponent_clauset(:,k),'r*-',...
            max_bouts_vector,wake_bout_exponent_MLE(:,k),'k*-');
        text(max_bouts_vector(end)+1,wake_bout_exponent_lin_reg(end,k),num2str(min_bouts_vector(k)));
        text(max_bouts_vector(end)+1,wake_bout_exponent_clauset(end,k),num2str(min_bouts_vector(k)));
        text(max_bouts_vector(end)+1,wake_bout_exponent_MLE(end,k),num2str(min_bouts_vector(k)));
    end
    title(['Wake bout exponents for ' sim_title]);
    legend('Linear reg','Clauset','MLE','Location','Best');
    hold off;
    
    if plot_full
        quotient = 1;
        [bsn,bsx] = hist(tempSb,round((max(tempSb)-min(tempSb))/quotient));
        [bwn,bwx] = hist(tempWb,round((max(tempWb)-min(tempWb))/quotient));
        
        figure
        subplot(2,2,1), loglog(bsx,bsn); title(sim_title);
        subplot(2,2,2), loglog(bwx,bwn);
        subplot(2,2,3), loglog(bsx,bsn); xlim([0.1 max(max_bouts_vector)]);
        subplot(2,2,4), loglog(bwx,bwn); xlim([0.1 max(max_bouts_vector)]);
    end
    
elseif detect_power_law_exponent_by_varying_min_bout_size
    which_experiment = 2;
    
    bout_definition = which_bout_definition;
    which_tpb_method = 'intervals are bouts';%'transition intervals are added to previous bout';%
    tic
    if ~use_current_workspace
        sleep_bouts = experiments(which_experiment).sleep_bouts;
        wake_bouts = experiments(which_experiment).wake_bouts;
        
        sim_title = experiments(which_experiment).info;
        A = experiments(which_experiment).A;
        WE = experiments(which_experiment).WE;
        SE = experiments(which_experiment).SE;
        nr_events = experiments(which_experiment).nr_events;
        t = experiments(which_experiment).t;
        
        [tempSb,tempWb,tempTb] = interval_to_bout(...
            experiments(which_experiment).duration_matrix,which_tpb_method);
    end
    
    quotient = 1;
    % Sleep bout exponent detection
    estimated_lower_min_sleep_bout = 0.1; %40;
    
    figure, hold on;
    title(['Sleep bout exponents for ' sim_title]);
    for j=0.1:0.05:1
        estimated_max_sleep_bout = mean(tempSb)*j;
        [min_bouts_vector, sleep_bout_exps] = detect_power_law_exp_network_case(...
            tempSb,quotient,...,
            estimated_lower_min_sleep_bout,estimated_max_sleep_bout);
        
        plot(min_bouts_vector,sleep_bout_exps);
        text(min_bouts_vector(end)+1,sleep_bout_exps(end),...
            num2str(round(estimated_max_sleep_bout)));
    end
    hold off;
    
    quotient = 4;
    [bsn,bsx] = hist(tempSb,round((max(tempSb)-min(tempSb))/quotient));
    [bwn,bwx] = hist(tempWb,round((max(tempWb)-min(tempWb))/quotient));
    
    figure
    subplot(2,2,1), loglog(bsx,bsn); title(sim_title);
    subplot(2,2,2), loglog(bwx,bwn);
    subplot(2,2,3), loglog(bsx,bsn);
    subplot(2,2,4), loglog(bwx,bwn);
elseif detect_power_law_exponent_automatically
    which_experiment = 3;
    which_tpb_method = 'transition intervals are added to previous bout';
    %'intervals are bouts';
    
    [tempSb,tempWb,tempTb] = interval_to_bout(...
        experiments(which_experiment).duration_matrix,which_tpb_method);
        
    quotient = 4;
    [bsn,bsx] = hist(tempSb,round((max(tempSb)-min(tempSb))/quotient));
    sleep_bout_title = ['Sleep bout distribution for ' ...
        experiments(which_experiment).info];
    
    figure,
    loglog(bsx,bsn); title(sleep_bout_title);
    

    tic
    [sleep_pl_lower_cutoff,sleep_pl_upper_cutoff,sleep_exponent] = ...
        detect_power_law_regime(tempSb,20:2:60,100:10:400,0.2,sleep_bout_title);
    toc
    
    % Plot fitted power-law on log-log bout distribution
    quotient1 = 4;
    figure,
    loglog(bsx,bsn); title(sleep_bout_title);
    
    hold on;
    x0 = sleep_pl_lower_cutoff; y0 = bsn(find(bsx>x0,1));
    x1 = sleep_pl_upper_cutoff; slope = sleep_exponent;
    plot([x0 x1],[y0,y0*(x1/x0)^(slope)],'ro-','Linewidth',2);
    hold off;
    legend('Pdf of sample data',...
        ['Fitted power-law in [' num2str(x0) ',' num2str(x1) ...
        '] w/ alpha=' num2str(slope)], ...
        'Location','Best');
    
end
clearvars -except experiments reg_bm_exps