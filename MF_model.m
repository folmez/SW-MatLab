function [t_out y_out] = MF_model(varargin)
% MF_MODEL
% Description of the system: \emph{states of neuron: $B$ (basal), $E$
% (excited), $I$ (inhibited) generate spikes at rate $\lambda_B, \lambda_E$
% and $\lambda_I$, respectively. Reception of excitatory spike $I \to B \to
% E$, reception of inhibitory spike $E \to B \to I$. Revert to basal state;
% $E\to B$ with rate $r_E$, $I \to B$ with rate $r_I$.}
%
%
% Let $w_B(t), w_I(t)$ and $w_E(t)$ be the fractions of wake neurons in
% basal, inhibited and excited states, respectively. Define $s_B(t),
% s_I(t)$ and $s_E(t)$ similarly.
%
% Define $\lambda_w(t):= w_E(t)\lambda_E +w_B(t)\lambda_B
% +w_I(t)\lambda_I$, the average spiking rate of neurons in the wake
% cluster. Define $\lambda_s(t):= s_E(t)\lambda_E +s_B(t)\lambda_B
% +s_I(t)\lambda_I$ similarly.
%
% Define $d_w$ and $d_s$ to be average degrees in wake and sleep clusters,
% respectively. These are counts of the average number of exicatory links a
% neuron has in each cluster. Average number of inhibitory links that a
% neuron has is defined to be $d_{ws}$. Notice that if exicatory links are
% assumed to be directed, then the averages $d_w$ and $d_s$ need only be
% halved to obtain average number of incoming excitatory links. For
% inhibitory links, we would need $d_{w\neg s}$ and $d_{s\neg w}$. I won't
% condider that case.

w_i_init = varargin{1};
w_e_init = varargin{2};
s_i_init = varargin{3};
s_e_init = varargin{4};
dW = varargin{5};
dS = varargin{6};
dIW = varargin{7};
dIS = varargin{8};
lambda_i = varargin{9};
lambda_b = varargin{10};
lambda_e = varargin{11};
ri = varargin{12};
rE = varargin{13};
t = varargin{14};
WE = varargin{15};
SE = varargin{16};
nW  = varargin{17}(1);
nS  = varargin{17}(2);
plot_simulation = varargin{18};
sim_title = varargin{19};
which_process = varargin{20};
% -------------------------------------------------------------------------
d_w = dW;
d_s = dS;

d_ws = dIW; % mean number of outgoing inhibitory degree of a wake node
d_sw = dIS; % mean number of outgoing inhibitory degree of a sleep node

t_end = 1e5;

y0 =[w_i_init; w_e_init; s_i_init; s_e_init];
tspan = [0 t_end];

[t_out, y_out] = ode45(@RHS,tspan,y0);

% Plot simulations and MF models on together
t_MF = find(t_out<max(t));

if plot_simulation
    figure
    plot(t, WE/nW, 'b', t, SE/nS, 'g', ...
    t_out(t_MF), y_out(t_MF,2), 'r', t_out(t_MF), y_out(t_MF,4), 'r');
    title(sim_title, 'FontSize', 15)
    h_xlabel = xlabel('Simulation time'); set(h_xlabel,'FontSize',15);
    h_ylabel = ylabel('Fraction of excited nodes'); set(h_ylabel,'FontSize',15);
    h_legend = legend('WE', 'SE','MF-Equilibriums','Location','Best');
    set(h_legend,'FontSize',15);
end
% -------------------------------------------------------------------------
    function dy = RHS(t,y)
        lambda_w = y(2)*lambda_e + (1-y(1)-y(2))*lambda_b + y(1)*lambda_i;
        lambda_s = y(4)*lambda_e + (1-y(3)-y(4))*lambda_b + y(3)*lambda_i;
        
        if strcmp(which_process,'three-state')
            dy1 = -ri*y(1) - (lambda_w)*d_w*y(1) + (lambda_s)*d_sw*(1-y(1)-y(2)); %WI
            dy2 = -rE*y(2) - (lambda_s)*d_sw*y(2) + (lambda_w)*d_w*(1-y(1)-y(2)); %WE
            dy3 = -ri*y(3) - (lambda_s)*d_s*y(3) + (lambda_w)*d_ws*(1-y(3)-y(4)); %SI
            dy4 = -rE*y(4) - (lambda_w)*d_ws*y(4) + (lambda_s)*d_s*(1-y(3)-y(4)); %SE
        elseif strcmp(which_process,'two-state')
            dy2 = -rE*y(2) - (lambda_s)*d_sw*y(2) + (lambda_w)*d_w*y(1); %WE
            dy4 = -rE*y(4) - (lambda_w)*d_ws*y(4) + (lambda_s)*d_s*y(3); %SE
            % Derivatives of WI(t) and SI(t) are negatives of WE(t) and SE(t)
            dy1 = -dy2;
            dy3 = -dy4;
        end
        
        dy = [dy1;dy2;dy3;dy4];
    end
% -------------------------------------------------------------------------
end