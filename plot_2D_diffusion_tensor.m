function plot_2D_diffusion_tensor(varargin)
% PLOT_2D_DIFFUSION_TENSOR plots eigenvalues of the diffusion tensor
% calculated for the stocashtic process (WE(t),SE(t))
%
%   e.g. plot_2D_diffusion_tensor(D, sim_title);
% ------------------------------------------------------------------------
D = varargin{1};    % Matrix of size 2 x 2 xn1+1 x n2+1
                    % (a matrix for each state point)
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
% In-script parameters
scaling=1;
% ------------------------------------------------------------------------
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
% Compute and plot diffusion tensor
figure, axis([0 n1 0 n2]);
title(sim_title, 'FontSize', 15);
h_xlabel = xlabel('SE'); set(h_xlabel,'FontSize',15);
h_ylabel = ylabel('WE'); set(h_ylabel,'FontSize',15);
hold on;
for i=1:k
    for j=1:k
        % Compute
        [eigV(:,:,i,j),eigD(:,:,i,j)] = eig(mean(mean( ...
            D(:, :, 1+(i-1)*n1/k:i*n1/k, ...
            1+(j-1)*n2/k:j*n2/k),3),4));
        vv = eigV(:,:,i,j); % eigenvectors: vv(:,1) and vv(:,2)
        dd = eigD(:,:,i,j); % eigenvalues: dd(1,1) and dd(2,2)
        % Plot
        plot([(i-0.5)*n1/k-scaling*dd(1,1)/2*vv(1,1),...
            (i-0.5)*n1/k+scaling*dd(1,1)/2*vv(1,1)],...
            [(j-0.5)*n2/k-scaling*dd(1,1)/2*vv(2,1),...
            (j-0.5)*n2/k+scaling*dd(1,1)/2*vv(2,1)]);
        plot([(i-0.5)*n1/k-scaling*dd(2,2)/2*vv(1,2),...
            (i-0.5)*n1/k+scaling*dd(2,2)/2*vv(1,2)],...
            [(j-0.5)*n2/k-scaling*dd(2,2)/2*vv(2,2),...
            (j-0.5)*n2/k+scaling*dd(2,2)/2*vv(2,2)]);
    end
end
h_legend = legend('Diffusion tensor','Location','Best');
set(h_legend,'FontSize',15);
axis([0 n1 0 n2]);
axis square;
hold off;
% ------------------------------------------------------------------------
end