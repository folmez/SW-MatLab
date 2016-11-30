function plot_activity_domains_on_surface_and_contour(varargin)
% PLOT_ACTIVITY_DOMAINS_ON_SURFACE_AND_CONTOUR(exp) plots approximated PDF
% of (WE,SE) as a surface and marks the activity domains on the surface. It
% also plots activity domains on a contour plot.
%
% Input: 'exp' must contain following cells:    A, 
%                                               wake_active_domain,
%                                               sleep_active_domain,
%                                               nW, nS or N and sim_title
% ------------------------------------------------------------------------
exp = varargin{1};
mark_domain_pts = 1;
want_surface = 1;
want_contour = 1;
plot_on_new_fig = 1;
fontsize = 15;
% ------------------------------------------------------------------------
i=2;
while i<=length(varargin),
    switch varargin{i},
        case 'exp',             exp = varargin{i+1};
        case 'mark_domain_pts', mark_domain_pts = varargin{i+1};
        case 'want_contour',    want_contour = varargin{i+1};
        case 'want_surface',    want_surface = varargin{i+1};
        case 'plot_on_new_fig', plot_on_new_fig = varargin{i+1};
        case 'FontSize',        fontsize = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
A = exp.WE_SE.density;
if mark_domain_pts
    % we need to add 1 in order to plot using plot3
    wad = exp.wake_active_domain + 1; 
    sad = exp.sleep_active_domain + 1;
end
if isfield(exp,'nW') && isfield(exp,'nS'),  nW = exp.nW;    nS = exp.nS;
elseif isfield(exp,'N')                     nW = exp.N;     nS = exp.N;
else                               [nW, nS] = size(A); nW=nW-1; nS=nS-1;
end
sim_title = exp.sim_title;

% Plot wake-active and sleep-active domains on surface
if want_surface
    if plot_on_new_fig, figure, end
    surf(A'), % transpose becase rows of A represent WE which needs to be the first-coordinate
    if mark_domain_pts
        hold,
        plot3(wad(:,1),wad(:,2),A(sub2ind([nW+1 nS+1],wad(:,1),wad(:,2))),'ko');
        % ,...    wpSE,wpWE,A(sub2ind([N+1 N+1],wpWE,wpSE)),'b*');
        plot3(sad(:,1),sad(:,2),A(sub2ind([nW+1 nS+1],sad(:,1),sad(:,2))),'kd');
        %, ...    spSE,spWE,A(sub2ind([N+1 N+1],spWE,spSE)),'b*');
    end
    axis tight;
    h_xlabel = xlabel('W_E'); set(h_xlabel, 'FontSize', fontsize);
    h_ylabel = ylabel('S_E'); set(h_ylabel, 'FontSize', fontsize);
    title(sim_title, 'FontSize', fontsize);
end
% Plot wake-active and sleep-active domains on contour
if want_contour
    if plot_on_new_fig, figure, end
    contour(A',100);
    if mark_domain_pts
        hold,
        plot(wad(:,1), wad(:,2), 'ko');
        plot(sad(:,1), sad(:,2), 'kd');
    end
    h_xlabel = xlabel('W_E'); set(h_xlabel, 'FontSize', fontsize);
    h_ylabel = ylabel('S_E'); set(h_ylabel, 'FontSize', fontsize);
    title(sim_title, 'FontSize', fontsize);
end
% ------------------------------------------------------------------------
end
