function varargout = do_SW_drift_diffusion_analysis(varargin)
%% Input arguments
nW = varargin{1}(1);
nS = varargin{1}(2);
sim_title = varargin{2};
t = varargin{3};
WE = varargin{4};
SE = varargin{5};
calculate_drift_diffusion = varargin{6};

plot_DDE_stuff = 0;

%% Model
if calculate_drift_diffusion
    [F,D] = calc_drift_diffusion_potential('2D', [nW nS], ...
        sim_title, t, WE, SE);
    [F_WE,D_WE] = calc_drift_diffusion_potential('1D', nW, ...
        sim_title, t, WE, 'WE', 'plot_DDE_stuff', plot_DDE_stuff);
    [F_SE,D_SE] = calc_drift_diffusion_potential('1D', nS, ...
        sim_title, t, SE, 'SE', 'plot_DDE_stuff', plot_DDE_stuff);
else
    F = zeros(nW+1,nS+1,2,1);
    D = zeros(2,2,nW+1,nS+1);
    F_WE = zeros(nW+1,1);
    D_WE = zeros(nW+1,1);
    F_SE = zeros(nS+1,1);
    D_SE = zeros(nS+1,1);
end

%% Output arguments
varargout{1} = F;
varargout{2} = D;
varargout{3} = F_WE;
varargout{4} = D_WE;
varargout{5} = F_SE;
varargout{6} = D_SE;

end

%% CALCULATE DRIFT, DIFFUSION AND POTENTIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, D] = calc_drift_diffusion_potential(varargin)
% CALCULATE_DRIFT_DIFFUSION_POTENTIAL computes a
% drift-diffusion-estimation to the 1D or 2D stochastic processes.
% e.g. for WE(t) or (WE(t),SE(t))

%% Input parameters
how_many_Ds = varargin{1}; % '1D' or '2D'
display_progress = 0;
plot_DDE_stuff = 1;

switch how_many_Ds
    case '1D'
        N = varargin{2};
        sim_title = varargin{3};
        t = varargin{4};
        rv = varargin{5}; % SE or WE
        data_title = varargin{6}; % 'SE' or 'WE'
    case '2D'
        n1 = varargin{2}(1);
        n2 = varargin{2}(2);
        sim_title = varargin{3};
        t = varargin{4};
        rv1 = varargin{5}; % SE
        rv2 = varargin{6}; % WE
    otherwise,
        display(varargin{1});
        error('Unexpected inputs!!!');
end

i=7;
while i<=length(varargin),
    switch varargin{i},
        case 'display_progress',        display_progress = varargin{i+1};
        case 'plot_DDE_stuff',          plot_DDE_stuff = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

%% Model
if strcmp(how_many_Ds, '2D')
    %% 2D DRIFT-DIFFUSION ESTIMATION
    % Drift vector and Diffusion tensor
    F = zeros(n1+1,n2+1,2,1); % row 1 and 2: drift in rv1 and rv2, resp.
    D = zeros(2,2,n1+1,n2+1);
    NaN_index = zeros((n1+1)*(n2+1),2); nan_count = 0;
    
    use_delta_t_from_data = 0;
    % 1: delta_t is computed from data
    % 0: delta_t is computed as mean(diff(t));
    
    delta_t = mean(diff(t(1:min(end,1e7))));
    
    tSIM = tic;
    for i=0:n1 % rv1
        for j=0:n2 % rv2
            % Calculating drift and diffusion
            x = find(rv1==i & rv2==j);
            
            if display_progress
                fprintf('(%2i,%2i):\t\t%i occurrences\n', i, j, length(x));
            end
            
            if ~isempty(x)
                if x(end)==length(rv1), x(end)=[]; end
                % Delta_rv: vector of displacements
                Delta_rv = [rv1(x+1)-rv1(x),rv2(x+1)-rv2(x)];
                delta_rv_tensor_delta_rv = zeros(2,2,length(x));
                for k=1:length(x)
                    % delta_rv: vector of displacements without the mean
                    %           displacement
                    
                    delta_rv = Delta_rv - ...
                        ones(length(x),1)*mean(Delta_rv,1);
                    delta_rv_tensor_delta_rv(:,:,k) = ...
                        delta_rv(k,:)'*delta_rv(k,:);
                end
                
                F(i+1,j+1,:,1) = mean(Delta_rv./delta_t,1);
                D(:,:,i+1,j+1) = 0.5 * ...
                    mean(delta_rv_tensor_delta_rv./delta_t,3);
            else
                nan_count = nan_count+1;
                NaN_index(nan_count,:) = [i+1 j+1];
            end
        end
    end
    fprintf('\n 2D DDE took %3.2f minutes\n\n', toc(tSIM)/60);
    
    % Plotting results
    if plot_DDE_stuff
        plot_2D_drift_vector(F, sim_title);
        plot_2D_diffusion_tensor(D, sim_title);
        plot_2D_peclet_number(F, D, sim_title);
    end
    
    fprintf('\n Drift-diffusion estimation took %3.2f minutes\n', ...
        toc(tSIM)/60);
    
elseif strcmp(how_many_Ds,'1D')
    %% 1D DRIFT-DIFFUSION ESTIMATION
    % Drift coefficient: F=F(WE+1)
    F = zeros(N+1,1);
    % Diffusion coefficient
    D = zeros(N+1,1);
    % Potential
    phi = zeros(N+1,1);
    peclet = zeros(N+1,1);
    
    use_delta_t_from_data = 0;
    % 1: delta_t is computed from data
    % 0: delta_t is computed as mean(diff(t));
    
    delta_t = mean(diff(t(1:min(end,1e6))));
    % delta_t = mean(diff(t));
    
    % Note that the frequency matrix is (N+1)x(N+1)
    tSIM = tic;
    for i=0:N
        % Calculating drift and diffusion
        x = find(rv==i);
        if ~isempty(x)
            if x(end)==length(rv), x(end)=[]; end
            
            Delta_rv = rv(x+1)-rv(x);
            if use_delta_t_from_data, delta_t = t(x+1)-t(x); end
            
            F(i+1) = mean(Delta_rv./delta_t);
            D(i+1) = 0.5 * mean((Delta_rv.^2)./delta_t);
            % [delta_rv(1:min(end,50)) temp(1:min(end,50))]
        end
        % Calculating potential
        if i>0
            X=0:i;
            phi(i+1) = (-1)*trapz(X,F(X+1)./D(X+1)) + log(D(i+1));
        end
        % Calculating Peclet number (Unitless number to compare drift and
        % diffusion)
        peclet(i+1) = (N) * F(i+1)/D(i+1);
    end
    fprintf('\n 1D DDE of %s took %3.2f minutes\n\n', ...
        data_title ,toc(tSIM)/60);
    
    % Plotting results
    [n,x] = hist(rv,N+1);
    [y1,max_ind_1] = max(n(1:N/2));
    [y2,max_ind_2] = max([zeros(1,N/2-1), n(N/2:N)]);
    plot_range = max_ind_1-1:max_ind_2+1;
    
    if plot_DDE_stuff
        figure,
        subplot(2,1,1); plot(x(plot_range),n(plot_range),'*');
        title(['Histogram of ' data_title ' from ' sim_title], ...
            'FontSize', 15);
        xlabel(data_title); xlim([min(plot_range) max(plot_range)]);
        
        subplot(2,1,2); plot(plot_range,peclet(plot_range));
        title(['Peclet number = N*(drift coeff)/(diffusion coeff) for ' ...
            data_title ' from ' sim_title], 'FontSize', 15);
        xlabel(data_title); axis([min(plot_range) max(plot_range) -5 5]);
        
        figure, subplot(2,1,1); plot(x,n,'*');
        title(['Histogram of ' data_title ' from ' sim_title], ...
            'FontSize', 15);
        xlabel(data_title);
        
        subplot(2,1,2); plot(0:N,phi);
        title(['Potential for ' data_title], 'FontSize', 15);
        xlabel(data_title);
        
        figure,
        subplot(2,1,1); plot(x,n,'*');
        title(['Histogram of ' data_title ' from ' sim_title], ...
            'FontSize', 15);
        xlabel(data_title);
        
        subplot(2,1,2); plot(0:N,F,0:N,D);
        title('Drift-Diffusion estimation', 'FontSize', 15)
        legend(['Drift coeff for ' data_title], ...
            ['Diffusion coeff for' data_title],'Location','Best');
        xlabel(data_title);
        hold on;
        plot([0 N],[0 0],'k.-.');
        hold off;
        
        X = 1:min(1e7,length(t)); % there is no need to plot all of t!
        figure,
        subplot(1,2,1), plot(X,t(X));
        xlabel('n'); ylabel('t(n)');
        title(['Mean delta t = mean(t(n+1)-t(n)) = ' ...
            num2str(delta_t)], 'FontSize', 15);
        
        [n_delta_t,x_delta_t] = hist(diff(t(X)),1000);
        
        subplot(1,2,2), plot(x_delta_t,log(n_delta_t),'b.-');
        hold on;
        title('Histogram of t(n+1)-t(n) = delta t', 'FontSize', 15);
        
        [n,x] = hist(exprnd(delta_t,length(X),1),1000);
        
        semilogy(x,log(n),'r.-');
        hold off;
        legend('Hist for simulated t(n+1)-t(n)', ...
            ['Hist for exp distributed data with mean' ...
            num2str(delta_t)], 'Location','Best');
        xlabel('delta t');
        ylabel('Natural logarithm of frequency');
    end
end

end