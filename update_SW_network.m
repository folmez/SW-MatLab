function [W, S, wait_time] = update_two_network_process(varargin)
%% Input parameters
i = varargin{1};            % Simulation event number
W = varargin{2};            % Wake network state vector
S = varargin{3};            % Sleep network state vector
WI = varargin{4};           % Number of (I)nhibited nodes in wake network
WE = varargin{5};           % Number of (E)xcited nodes in wake network
SI = varargin{6};           % Number of (I)nhibited nodes in sleep network
SE = varargin{7};           % Number of (E)xcited nodes in sleep network
WW = varargin{8};           % Wake network adjacency matrix
SS = varargin{9};           % Sleep network adjacency matrix
WS = varargin{10};          % Wake to Sleep inhibitory connections
SW = varargin{11};          % Sleep to Wake inhibitory connections
nW = varargin{12}(1);       % Size of nodes in each network
nS = varargin{12}(2);   
lambda_i = varargin{13};    % Firing rate of an (I)nhibited node
lambda_b = varargin{14};    % Firing rate of an (B)asal node
lambda_e = varargin{15};    % Firing rate of an (E)xcited node
ri = varargin{16};          % Relaxation rate of an (I)nhibited node
rE = varargin{17};          % Relaxation rate of an (E)xcited node
use_gillespie_algorithm = varargin{18};
display_sim_summary = varargin{19};
which_process = varargin{20};   % 'two-state' or 'three-state'
%------------------------------------------------------------------
display_details = 0;
if strcmp(which_process,'three-state')
    bump_step_size = 1;
elseif strcmp(which_process,'two-state')
    bump_step_size = 2;
end
%------------------------------------------------------------------
%% Model
if use_gillespie_algorithm
    total_rate = WI(i-1)*ri + WE(i-1)*rE + SI(i-1)*ri + SE(i-1)*rE +...
        WI(i-1)*lambda_i + (nW-WI(i-1)-WE(i-1))*lambda_b + WE(i-1)*lambda_e +...
        SI(i-1)*lambda_i + (nS-SI(i-1)-SE(i-1))*lambda_b + SE(i-1)*lambda_e;
    % Total_rate
    wait_time = exprnd(1/total_rate);
    prob_events = [(W==-1)*lambda_i + (W==0)*lambda_b + (W==1)*lambda_e;...
        (S==-1)*lambda_i + (S==0)*lambda_b + (S==1)*lambda_e;...
        (W==-1)*ri + (W==1)*rE; (S==-1)*ri + (S==1)*rE]/total_rate;
    cum_prob_events = cumsum(prob_events);
    cum_prob_events(prob_events==0) = 0;
    % cum_prob_events
    toss = rand;
    % [cum_prob_events cum_prob_events<=toss]
    [~,I] = min(cum_prob_events<=toss);
else
    W_excitation_times = (W==-1).*exprnd(1/lambda_i,nW,1)+...
        (W==0).*exprnd(1/lambda_b,nW,1)+(W==1).*exprnd(1/lambda_e,nW,1);
    W_relaxation_times = (W==-1).*exprnd(1/ri,nW,1)+...
        (W==1).*exprnd(1/rE,nW,1); W_relaxation_times(~W_relaxation_times)=inf;
    S_relaxation_times = (S==-1).*exprnd(1/ri,nS,1)+...
        (S==1).*exprnd(1/rE,nS,1); S_relaxation_times(~S_relaxation_times)=inf;
    S_excitation_times = (S==-1).*exprnd(1/lambda_i,nS,1)+...
        (S==0).*exprnd(1/lambda_b,nS,1)+(S==1).*exprnd(1/lambda_e,nS,1);
    
    [wait_time,I] = min([W_excitation_times;S_excitation_times;...
        W_relaxation_times;S_relaxation_times]);
end
%------------------------------------------------------------------
%% Update process
if I > 2*nW + nS    % sleep node relaxes
    if display_sim_summary
        if S(I-(2*nW+nS))==-1
            fprintf('%i\t: I-sleep node #%i relaxed\n', i, I-(2*nW+nS));
        elseif S(I-(2*nW+nS))==1
            fprintf('%i\t: E-sleep node #%i relaxed\n', i, I-(2*nW+nS));
        end
    end
    S(I-(2*nW+nS)) = 0;
elseif I > nW + nS  % wake node relaxes
    if display_sim_summary
        if W(I-(nW+nS))==-1
            fprintf('%i\t: I-wake node #%i relaxed\n', i, I-(nW+nS));
        elseif W(I-(nW+nS))==1
            fprintf('%i\t: E-wake node #%i relaxed\n', i, I-(nW+nS));
        end
    end
    W(I-(nW+nS)) = 0;
elseif I > nW       % sleep node fires
    if display_sim_summary
        if S(I-nW)==-1
            fprintf('%i\t: I-sleep node #%i FIRED\n', i, I-nW);
        elseif S(I-nW)==0
            fprintf('%i\t: B-sleep node #%i FIRED\n', i, I-nW);
        elseif S(I-nW)==1
            fprintf('%i\t: E-sleep node #%i FIRED\n', i, I-nW);
        end
    end
    S = S + bump_step_size * SS(I-nW,:)';
    S(S>1)=1;
    W = W - bump_step_size * SW(I-nW,:)';
    W(W<-1)=-1;
else                % wake node fires
    if display_sim_summary
        if W(I)==-1
            fprintf('%i\t: I-wake node #%i FIRED\n', i, I);
        elseif W(I)==0
            fprintf('%i\t: B-wake node #%i FIRED\n', i, I);
        elseif W(I)==1
            fprintf('%i\t: E-wake node #%i FIRED\n', i, I);
        end
    end
    W = W + bump_step_size * WW(I,:)';
    W(W>1)=1;
    S = S - bump_step_size * WS(I,:)';
    S(S<-1)=-1;
end
%% Display results
if display_details
    display(W);
    display(WW);
    display(WS);
    display(S);
    display(SS);
    display(SW);
    pause;
end
%------------------------------------------------------------------
end