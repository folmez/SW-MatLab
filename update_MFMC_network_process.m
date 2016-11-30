function [WE_next, WI_next, SE_next, SI_next, t_next] = ...
    update_MFMC_network_process(varargin)
% UPDATE_MFMC_NETWORK_PROCESS updates a given mean field Markov chain
% process on two networks.
%
i = varargin{1};            % Simulation event number
WE_prev = varargin{2};      % Previous count of Wake-Excited nodes
WI_prev = varargin{3};      % Previous count of Wake-Inhibited nodes
SE_prev = varargin{4};      % Previous count of Sleep-Excited nodes
SI_prev = varargin{5};      % Previous count of Sleep-Inhibited nodes
t_prev = varargin{6};      % Previous total time

EdegW = varargin{7};        % Excitatory degree dist of wake cluster
EdegS = varargin{8};        % Excitatory degree dist of sleep cluster
IdegW = varargin{9};        % Inhibitory degree dist of wake cluster
IdegS = varargin{10};        % Inhibitory degree dist of sleep cluster

N = varargin{11};           % Size of nodes in each network

lambda_i = varargin{12};    % Firing rate of an (I)nhibited node
lambda_b = varargin{13};    % Firing rate of an (B)asal node
lambda_e = varargin{14};    % Firing rate of an (E)xcited node
ri = varargin{15};          % Relaxation rate of an (I)nhibited node
rE = varargin{16};          % Relaxation rate of an (E)xcited node

display_sim_summary = varargin{17};
%------------------------------------------------------------------
% [wait_time,I] = min([...
%     exprnd(1/(WE_prev*lambda_e + (N-WE_prev-WI_prev)*lambda_b + WI_prev*lambda_i)),...
%     exprnd(1/(SE_prev*lambda_e + (N-SE_prev-SI_prev)*lambda_b + SI_prev*lambda_i)),...
%     exprnd(1/(WE_prev*rE + WI_prev*ri)),...
%     exprnd(1/(SE_prev*rE + SI_prev*ri))]);

WB_prev = N-WE_prev-WI_prev;
SB_prev = N-SE_prev-SI_prev;

% Gillespie algorithm
total_rate = WI_prev*ri + WE_prev*rE + SI_prev*ri + SE_prev*rE +...
    WI_prev*lambda_i + WB_prev*lambda_b + WE_prev*lambda_e +...
    SI_prev*lambda_i + SB_prev*lambda_b + SE_prev*lambda_e;
prob_events = [...
    WI_prev*lambda_i; WB_prev*lambda_b; WE_prev*lambda_e; ...
    SI_prev*lambda_i; SB_prev*lambda_b; SE_prev*lambda_e; ...
    WI_prev*ri; WE_prev*rE; ...
    SI_prev*ri; SE_prev*rE] / total_rate;
% WI fires, WB fires, WE fires, SI fires, SB fires, SE fires,
% WI relaxes, WE relaxes, SI relaxes, SE relaxes

cum_prob_events = cumsum(prob_events);
cum_prob_events(prob_events==0) = 0;
% cum_prob_events
toss = rand;
% [cum_prob_events cum_prob_events<=toss]
[~,I] = min(cum_prob_events<=toss);
% Wait time for next event
delta_t = exprnd(1/total_rate);

%------------------------------------------------------------------
% Update process
node_number = randi(N,1);
fWnid = EdegW(node_number); % firing Wake node in-degree
fWnod = IdegW(node_number); % firing Wake node out-degree
fSnid = EdegS(node_number); % firing Sleep node in-degree
fSnod = IdegS(node_number); % firing Sleep node out-degree

if SI_prev<0 || SE_prev<0 || WI_prev<0 ||WE_prev<0 || ...
        N-SI_prev-SE_prev<0 || N-WI_prev-WE_prev<0 
    display('darn!');
end


delta_WE = 0;
delta_WI = 0;
delta_SI = 0;
delta_SE = 0;

if I==1 % I-wake node fires
    delta_WE = min(binornd(fWnid,WB_prev/N), WB_prev);      % WB up to WE
    delta_WI = (-1)*binornd(fWnid-delta_WE,(WI_prev-1)/N);  % WI up to WB
    delta_SI = min(binornd(fWnod,SB_prev/N), SB_prev);      % SB down to SI
    delta_SE = (-1)*binornd(fWnod-delta_SI,SE_prev/N);      % SE down to SB 
elseif I==2 % B-wake node fires
    delta_WE = min(binornd(fWnid,(WB_prev-1)/N), WB_prev-1);% WB up to WE
    delta_WI = (-1)*binornd(fWnid-delta_WE,WI_prev/N);      % WI up to WB   
    delta_SI = min(binornd(fWnod,SB_prev/N), SB_prev);      % SB down to SI
    delta_SE = (-1)*binornd(fWnod-delta_SI,SE_prev/N);      % SE downto SB
elseif I==3 % E-wake node fires
    delta_WE = min(binornd(fWnid,WB_prev/N), WB_prev);      % WB up to WE
    delta_WI = (-1)*binornd(fWnid-delta_WE,WI_prev/N);      % WI up to WB
    delta_SI = min(binornd(fWnod,SB_prev/N), SB_prev);      % SB down to SI
    delta_SE = (-1)*binornd(fWnod-delta_SI,SE_prev/N);      % SE down to SB
elseif I==4 % I-sleep node fires
    delta_WI = min(binornd(fSnod,WB_prev/N), WB_prev);      % WB down to WI
    delta_WE = (-1)*binornd(fSnod+delta_WI,WE_prev/N);      % WE down to WB
    delta_SE = min(binornd(fSnid,SB_prev/N), SB_prev);      % SB up to SE
    delta_SI = (-1)*binornd(fSnid-delta_SE,(SI_prev-1)/N);  % SI up to SB
elseif I==5 % B-sleep node fires
    delta_WI = min(binornd(fSnod,WB_prev/N), WB_prev);      % WB down to WI
    delta_WE = (-1)*binornd(fSnod+delta_WI,WE_prev/N);      % WE down to WB
    delta_SE = min(binornd(fSnid,(SB_prev-1)/N), SB_prev-1);% SB up to SE
    delta_SI = (-1)*binornd(fSnid-delta_SE,SI_prev/N);      % SI up to SB
elseif I==6 % E-sleep node fires
    delta_WI = min(binornd(fSnod,WB_prev/N), WB_prev);      % WB down to WI
    delta_WE = (-1)*binornd(fSnod+delta_WI,WE_prev/N);      % WE down to WB
    delta_SE = min(binornd(fSnid,SB_prev/N), SB_prev);      % SB up to SE
    delta_SI = (-1)*binornd(fSnid-delta_SE,SI_prev/N);      % SI up to SB
elseif I==7 % I-wake node relaxes
    delta_WI = -1;
elseif I==8 % E-wake node relaxes
    delta_WE = -1;
elseif I==9 % I-sleep node relaxes
    delta_SI = -1;
elseif I==10 % E-sleep node relaxes
    delta_SE = -1;
else
    error('A problem occurred while updating MFMC two network process!!!');
end

WE_next = median([0; N; WE_prev + delta_WE]);
SE_next = median([0; N; SE_prev + delta_SE]);
WI_next = median([0; N; WI_prev + delta_WI]);
SI_next = median([0; N; SI_prev + delta_SI]);
if  N-SI_next-SE_next<0 || N-WI_next-WE_next<0;
    display('darn!');
end
t_next = t_prev + delta_t;
%------------------------------------------------------------------
% Display update summary
if display_sim_summary
    if I==1,        fprintf('%i\t: I-wake node fired',i);
    elseif I==2,    fprintf('%i\t: B-wake node fired',i);
    elseif I==3,    fprintf('%i\t: E-wake node fired',i);
    elseif I==4,    fprintf('%i\t: I-sleep node fired',i);
    elseif I==5,    fprintf('%i\t: B-sleep node fired',i);
    elseif I==6,    fprintf('%i\t: E-sleep node fired',i);
    elseif I==7,    fprintf('%i\t: I-wake node relaxed',i);
    elseif I==8,    fprintf('%i\t: E-wake node relaxed',i);
    elseif I==9,    fprintf('%i\t: I-sleep node relaxed',i);
    elseif I==10,   fprintf('%i\t: E-sleep node relaxed',i);
    end
    fprintf('\t%i\t%i\t%i\t%i\t%i\t%i\t\n',WE_next, N-WE_next-WI_next, ...
        WI_next, SE_next, N-SE_next-SI_next, SI_next);
end
%------------------------------------------------------------------
end