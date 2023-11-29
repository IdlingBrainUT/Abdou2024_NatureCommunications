%% approximative_temporal_binning_code
function [approx_t, approx_val, approx_binned_data] = fxn_approx_temporal_binning(temp_time, temp_value, binning_unit)
%% for debug
% temp_time         =  results_yd1.total_Stage_time_val(:,1);
% temp_value        =  results_yd1.total_Stage_time_val(:,2);
% binning_unit = 0.2 ; % second 
%% time info
approx_temp_time_end = ceil(temp_time(end));
approx_t = [binning_unit: binning_unit: approx_temp_time_end];
%%
search_cell_temp = [];
% for i = 1 % for debug   
for i = 1:size(approx_t,2)
    if i <= size(approx_t,2)-1
search_onset = approx_t(i);
search_end   = approx_t(i+1);
search_cell_temp  = find(search_onset <= temp_time & temp_time < search_end); %  < search_end)
        if sum(search_cell_temp) == 0
%             disp('stop cal')
        else
        approx_val(i) = max(temp_value(search_cell_temp)); 
        end
    elseif i == size(approx_t,2)
search_onset = approx_t(i-1);
search_end   = approx_t(i);
search_cell_temp  = find(search_onset <= temp_time & temp_time < search_end); %  < search_end)
        if sum(search_cell_temp) == 0
%             disp('stop cal')
        else
        approx_val(i) = max(temp_value(search_cell_temp)); 
        end
    end
end
%%
if numel(approx_t) ~= numel(approx_val)
    approx_t = approx_t(1:numel(approx_val));
disp('aprrox_val is truncated');
else
disp('normally processed.')
end
%%
approx_t_tp   = approx_t';
approx_val_tp = approx_val';
approx_binned_data = [approx_t_tp, approx_val_tp];
disp('Done approximate temporal binning')
%% figure
% figure; plot(approx_binned_data(:,1), approx_binned_data(:,2));

%% sol time stamp
% temp_sol                  = total_SOL_time_val(:,1:2);
% temp_sol_timestamp_frame  = find(temp_sol(:,2)==1); % frame num
% temp_sol_timestamp_time   = temp_sol(temp_sol_timestamp_frame,1); % in second
%% alighment between sol and spectrogram stamps
% % knnsearch(X,Y); % X: data set, Y: target, Nearest neighbor search
% % Y 内の各クエリ点に対する最近傍を X 内で探索し、列ベクトル Idx にインデックスを返します。Idx の行数は Y と同じです。
% for i = 1:size(temp_sol_timestamp_time,1)
% temp_MPPCA_t_knn_frame(i) = knnsearch(MPPCA_t', temp_sol_timestamp_time(i));
% end
% temp_MPPCA_timestamp_time = MPPCA_t(temp_MPPCA_t_knn_frame);
% temp_timestamp_master = [temp_sol_timestamp_frame, temp_sol_timestamp_time, temp_MPPCA_t_knn_frame', temp_MPPCA_timestamp_time'];
%%
end


