function plot_drift_diffusion_results(varargin)
% PLOT_DRIFT_DIFFUSION_RESULTS(experiment) plots already calculated
% drift-diffusion results
% ------------------------------------------------------------------------
i=1;
while i<=length(varargin),
    switch varargin{i},
        case 'experiment'
            experiment = varargin{i+1};
        case 'what_DDE'
            what_DDE = varargin{i+1}; % '2D' or 'WE' or 'SE'
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
switch what_DDE
    case '2D'
        sim_title = experiment.sim_title;
        F = experiment.F;
        D = experiment.D;
        if isfield(experiment,'nW')
            n1 = experiment.nW;
            n2 = experiment.nS;
        elseif isfield(experiment,'N')
            n1 = experiment.N;
            n2 = experiment.N;
        end
        
        U = F(:,:,1,1);
        V = F(:,:,2,1);
        
        k = 20;
        % Plot drift coefficient
        U_new=zeros(k,k);
        V_new=zeros(k,k);
        for i=1:k
            for j=1:k
                U_new(i,j)=sum(sum(U(1+(i-1)*n1/k:i*n1/k,1+(j-1)*n2/k:j*n2/k)));
                V_new(i,j)=sum(sum(V(1+(i-1)*n1/k:i*n1/k,1+(j-1)*n2/k:j*n2/k)));
            end
        end
        [X,Y] = meshgrid(n1/k:n1/k:n1,n2/k:n2/k:n2);
        figure,
        quiver(X,Y,U_new,V_new);
        title(['Drift vector field (' sim_title ')']);

        % -----------------------------------------
        % Plot average diffusion tensor
        eigV = zeros(2,2,k,k); % 4D - eigenvector array
        eigD = zeros(2,2,k,k); % 4D - eigenvalue array
        
        scaling = 1;
        figure, hold on;
        for i=1:k
            for j=1:k
                [eigV(:,:,i,j),eigD(:,:,i,j)] = ...
                    eig(mean(mean(D(:, :, 1+(i-1)*n1/k:i*n1/k, ...
                    1+(j-1)*n2/k:j*n2/k),3),4));
                
                vv = eigV(:,:,i,j); % eigenvectors: vv(:,1) and vv(:,2)
                dd = eigD(:,:,i,j); % eigenvalues: dd(1,1) and dd(2,2)
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
        axis([0 n1 0 n2]);
        title(['Diffusion tensors for (' sim_title ')']);
        xlabel('SE'); ylabel('WE');
        hold off;
        % -----------------------------------------
        % Plot peclet number contour map for 2D system
        
        % (1/2) Drift coefficient
        U_new = zeros(n1+1,n2+1);
        V_new = zeros(n1+1,n2+1);
        for i=1:n1+1
            for j=1:n2+1
                U_new(i,j)=sum(sum(U(1+(i-1):i,1+(j-1):j)));
                V_new(i,j)=sum(sum(V(1+(i-1):i,1+(j-1):j)));
            end
        end
        
        % (2/2) Diffusion tensor eigenvalues
        eigD = zeros(2,2,n1+1,n2+1); % 4D - eigenvalue array
        for i=1:n1+1
            for j=1:n2+1
                [~,eigD(:,:,i,j)] = eig(D(:,:,i,j));
                dd = eigD(:,:,i,j); % eigenvalues: dd(1,1) and dd(2,2)
            end
        end
        
        % Plot Peclet number for 2D system
        peclet = zeros(n1+1,n2+1);
        figure,
        for i=1:n1+1
            for j=1:n2+1
                dd = eigD(:,:,j,i); % eigenvalues: dd(1,1) and dd(2,2)
                peclet(i,j) = min(n1,n2)*... % length scale
                    norm([U_new(i,j) V_new(i,j)])/max(dd(1,1),dd(2,2));
            end
        end
        peclet_truncated = peclet;
        peclet_truncated(peclet>50) = 50;
        [c, h] = contour(peclet_truncated);
        clabel(c, h);
        title(['Peclet number for 2D drift-diffusion estimation of' sim_title]);
        xlabel('SE'); ylabel('WE');        
end
% ------------------------------------------------------------------------
end
