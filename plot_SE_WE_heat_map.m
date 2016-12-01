function plot_SE_WE_heat_map(varargin)
%% Input parameters
A = varargin{1};
sim_title = varargin{2};
plot_heat_map = varargin{3};
pause_after_plot = 0;

i = 4;
while i<=length(varargin),
    switch varargin{i},
        case 'pause_after_plot'
            pause_after_plot = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

%% Plot heat map
if plot_heat_map
    figure,
    surf(A);
    axis tight;
    xlabel('# of excited SLEEP nodes');
    ylabel('# of excited WAKE nodes');
    if pause_after_plot
        title('Press a key to continue if bi-stable...');
        drawnow;
        pause;
    else     
        title(sim_title);
    end
end

end