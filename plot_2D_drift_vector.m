function plot_2D_drift_vector(varargin)
% PLOT_2D_DRIFT_VECTOR plots drift vector calculated for the stocashtic
% process (WE(t),SE(t))
%
%   e.g. plot_2D_drift_vector(F, sim_title);
% ------------------------------------------------------------------------
F = varargin{1};    % Matrix of size n1+1 x n2+1 x 2 x 1
                    % (a vector for each state point)
sim_title = varargin{2};
nr_pts_to_show = 20;
% ------------------------------------------------------------------------
i=3;
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
U = F(:,:,1,1); V = F(:,:,2,1);
temp = size(F);
n1 = temp(1)-1;
n2 = temp(2)-1;
if mod(n1,nr_pts_to_show)~=0 || mod(n2,nr_pts_to_show)~=0
    error('Number of points desired is not a divisor of graph sizes');
end
k = nr_pts_to_show;
% ------------------------------------------------------------------------
drift_norm = zeros(n1+1, n2+1);
% Compute norm of drift
for i=1:n1+1
    for j=1:n2+1
        drift_norm(i,j) = norm([U(i,j) V(i,j)]);
    end
end

% Plot norm of drift
figure, surface(drift_norm);
axis tight;
colorbar('Location','EastOutSide');

% Compute drift vectors
U_new=zeros(k,k); V_new=zeros(k,k);
for i=1:k
    for j=1:k
        U_new(i,j)=sum(sum(U(1+(i-1)*n1/k:i*n1/k,1+(j-1)*n2/k:j*n2/k)));
        V_new(i,j)=sum(sum(V(1+(i-1)*n1/k:i*n1/k,1+(j-1)*n2/k:j*n2/k)));
    end
end

% Plot drift vectors using arrows
[X,Y] = meshgrid(n1/k:n1/k:n1,n2/k:n2/k:n2);
figure,
quiver(X,Y,U_new,V_new);
axis tight;
title(sim_title, 'FontSize', 15);
h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
h_legend = legend('Drift vector field','Location','Best');
set(h_legend,'FontSize',15);
% ------------------------------------------------------------------------
end