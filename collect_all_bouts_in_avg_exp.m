function avg_exp = collect_all_bouts_in_avg_exp(varargin)
exps = varargin{1};
avg_exp = varargin{2};
nr_total_bouts = 1e6;
% ------------------------------------------------------------------------
ALL_bouts = zeros(nr_total_bouts, 6); 
ALL_bouts_size = 0;
for i = 1:length(exps)
    % All bouts
    end_idx = min([length(exps(i).bouts.wake.no_trans.data) ...
        length(exps(i).bouts.sleep.no_trans.data)]);
    bouts_matrix = [exps(i).bouts.wake.no_trans.data(1:end_idx), ...
        exps(i).bouts.wake.prev_trans.data(1:end_idx), ...
        exps(i).bouts.wake.half_trans.data(1:end_idx), ...
        exps(i).bouts.sleep.no_trans.data(1:end_idx), ...
        exps(i).bouts.sleep.prev_trans.data(1:end_idx), ...
        exps(i).bouts.sleep.half_trans.data(1:end_idx)];
    nr_bouts = length(bouts_matrix);
    if ALL_bouts_size + nr_bouts <= length(ALL_bouts)
        ALL_bouts(ALL_bouts_size+1:ALL_bouts_size+nr_bouts, :) = bouts_matrix;
        ALL_bouts_size = ALL_bouts_size + nr_bouts;
    end
end
ALL_bouts(ALL_bouts_size+1:end,:)=[];
% ------------------------------------------------------------------------
ALL_WS_trans_ints = zeros(nr_total_bouts, 1);
ALL_WS_trans_ints_size = 0;
ALL_SW_trans_ints = zeros(nr_total_bouts, 1);
ALL_SW_trans_ints_size = 0;
for i = 1:length(exps)
    % All transition intervals
    [~,~,~, WSti, SWti] = interval_to_bout(exps(i).duration_matrix, ...
        'intervals are bouts');
    nr_WSti = length(WSti);
    nr_SWti = length(SWti);
    if ALL_WS_trans_ints_size + nr_WSti <= nr_total_bouts
        ALL_WS_trans_ints(ALL_WS_trans_ints_size+1:...
            ALL_WS_trans_ints_size+nr_WSti, :)= ...
            WSti;
        ALL_WS_trans_ints_size = ALL_WS_trans_ints_size + nr_WSti;
    end
    if ALL_SW_trans_ints_size + nr_SWti <= nr_total_bouts
        ALL_SW_trans_ints(ALL_SW_trans_ints_size+1:...
            ALL_SW_trans_ints_size+nr_SWti, :)= ...
            SWti;
        ALL_SW_trans_ints_size = ALL_SW_trans_ints_size + nr_SWti;
    end
end
ALL_WS_trans_ints(ALL_WS_trans_ints_size+1:end,:)=[];
ALL_SW_trans_ints(ALL_SW_trans_ints_size+1:end,:)=[];
% ------------------------------------------------------------------------
% Update avg_exp
avg_exp.bouts.wake.no_trans.data = ALL_bouts(:,1);
avg_exp.bouts.wake.prev_trans.data = ALL_bouts(:,2);
avg_exp.bouts.wake.half_trans.data = ALL_bouts(:,3);
avg_exp.bouts.sleep.no_trans.data = ALL_bouts(:,4);
avg_exp.bouts.sleep.prev_trans.data = ALL_bouts(:,5);
avg_exp.bouts.sleep.half_trans.data = ALL_bouts(:,6);

avg_exp.intervals.transition.WS = ALL_WS_trans_ints;
avg_exp.intervals.transition.SW = ALL_SW_trans_ints;
% ------------------------------------------------------------------------
end