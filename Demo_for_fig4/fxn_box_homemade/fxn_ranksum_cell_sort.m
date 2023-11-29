%% NREM REM phase judgement
% non-parametric analysis using Wilcoxon rank sum test
% [p,h,stats] = ranksum(mileage(:,1),mileage(:,2))
function [result_occ_ranksum, result_sort_cell_id, sorted_ca_data] = fxn_ranksum_cell_sort(data_temp, all_time_stamp, time_stamp, fold_thresold) % 
%% for debug,
% data_temp =  data_temp;
% time_stamp = time_stamp;
%% Calculate sorting frames
fig_clim_range = [0 3];
upsampling_factor = 1;
control_frame =  (1:size(data_temp,1))'; % copy
control_frame(time_stamp) = [];
test_frame = time_stamp;
%%
% for sorted fig
occ_sorted  = data_temp([control_frame;test_frame],:); % NREM REM sort
occ_control = data_temp(control_frame,:);
occ_test    = data_temp(time_stamp,:);
%% figure data
occ_temp_tp   = data_temp';
occ_sorted_tp = occ_sorted';
%%
results_rank = {}; results_rank{1,2} = 'Occ ctrl'; results_rank{1,3} = 'Occ test';
results_rank{1,4} = 'ctrl mean'; results_rank{1,5} = 'test mean'; results_rank{1,6} = 'Dominance';
results_rank{1,7} = 'Ranksum p'; results_rank{1,8} = 'ranksum h'; results_rank{1,9} = 'ranksum stats'; 
results_rank{1,10} = 'fold'; results_rank{1,11} = 'test sort'; 
sort_cell_id_temp = [];

for i = 1:size(occ_control,2)
    pattern_name = ['cell ID #', num2str(i)];
    results_rank{1+i,1} = pattern_name;
    results_rank{1+i,2} = occ_control(:,i);
    results_rank{1+i,3} = occ_test(:,i);
    results_rank{1+i,4} = mean(results_rank{1+i,2},1);
    results_rank{1+i,5} = mean(results_rank{1+i,3},1);
    
    [results_rank{1+i,7}, results_rank{1+i,8}, results_rank{1+i,9}] = ranksum(results_rank{1+i,2},results_rank{1+i,3});
    
    if     (results_rank{1+i,5} > results_rank{1+i,4}) && (results_rank{1+i,8} == 1)
        results_rank{1+i,6} = 'test';
    elseif (results_rank{1+i,5} < results_rank{1+i,4}) && (results_rank{1+i,8} == 1)
        results_rank{1+i,6} = 'control';
    else
        results_rank{1+i,6} = 'non sig';
    end
    
    results_rank{1+i,10} =  results_rank{1+i,5}/results_rank{1+i,4};
    results_rank{1+i,11} =  (results_rank{1+i,10} > fold_thresold) && (results_rank{1+i,8} == 1);
    
    sort_cell_id_temp =  cat(1, sort_cell_id_temp, double(results_rank{1+i,11})); 
end
%%
result_sort_cell_id = find(sort_cell_id_temp == 1);
result_occ_ranksum = results_rank;
%%

ID_sig_cell = result_sort_cell_id;
ID_others   = 1:size(data_temp,2);
ID_others(ID_sig_cell) = [];

sorted_ca_data  = data_temp(:,[ID_sig_cell; ID_others']); 
%%
figure
subplot(7,1,1:2); imagesc([1:size(occ_temp_tp,2)]./upsampling_factor, 1:size(occ_temp_tp,1), occ_temp_tp); 
clim(fig_clim_range); title('raw data')

subplot(7,1,3:4); imagesc([1:size(occ_sorted_tp,2)]./upsampling_factor, 1:size(occ_sorted_tp,1), occ_sorted_tp); 
clim(fig_clim_range); title('temporal sorting')

subplot(715); plot(all_time_stamp,'k');title('time stamp')
subplot(7,1,6:7); imagesc([1:size(occ_temp_tp,2)]./upsampling_factor, 1:size(occ_temp_tp,1), sorted_ca_data'); 
clim(fig_clim_range); title('cell ID sorted')
%% disp
disp(results_rank(:,[1 4 5 6 10 11]));
%%
end