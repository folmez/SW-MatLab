function plot_2D_peclet_number(varargin)
% PLOT_2D_PECLET_NUMBER plots a contour map of Peclet number calculated for
% the stocashtic process (WE(t),SE(t)) at every state points
%
%   e.g. plot_2D_peclet_number(F, D, sim_title);
% ------------------------------------------------------------------------
F = varargin{1};    % Matrix of size n1+1 x n2+1 x 2 x 1
                    % (a vector for each state point)
D = varargin{2};    % Matrix of size 2 x 2 xn1+1 x n2+1
                    % (a matrix for each state point)
sim_title = varargin{3};
nr_pts_to_show = 20;
% ------------------------------------------------------------------------
i=4;
while i<=length(varargin),
    switch varargin{i},
        case 'nr_pts_to_show',      nr_pts_to_show = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
% In-script parameters
scaling=1;
% ------------------------------------------------------------------------
U = F(:,:,1,1); V = F(:,:,2,1);
temp = size(D);
n1 = temp(3)-1;
n2 = temp(4)-1;
if mod(n1,nr_pts_to_show)~=0 || mod(n2,nr_pts_to_show)~=0
    error('Number of points desired is not a divisor of graph sizes');
end
k = nr_pts_to_show;
% ------------------------------------------------------------------------
eigV = zeros(2,2,n1+1,n2+1); % 4D - eigenvector array
eigD = zeros(2,2,n1+1,n2+1); % 4D - eigenvalue array
% Compute and plot Peclet number for 2D system
peclet = zeros(n1+1,n2+1);
for i=1:n1+1
    for j=1:n2+1
        % Eigenvalues: dd(1,1) and dd(2,2)
        [~ , dd] = eig(D(:,:,i,j));
        peclet(i,j) = max([n1 n2]) * ...
            norm([U(i,j) V(i,j)]) / max(dd(1,1),dd(2,2));
    end
end
figure,
surface(peclet);
title(sim_title, 'FontSize', 15);
h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
h_legend = legend('Peclet number','Location','Best');
set(h_legend,'FontSize',15);
axis tight;

peclet_truncated = peclet;
peclet_truncated(peclet>50) = 50;
figure,
[c, h] = contour(peclet_truncated);
clabel(c, h);
title(sim_title, 'FontSize', 15);
h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
h_legend = legend('Peclet number','Location','Best');
set(h_legend,'FontSize',15);
% ------------------------------------------------------------------------
end