%% MPPCA cell-cell pairwise interation 
function [result_multi] = fxn_MPPCA_cal_connectivity_multi_core(result_data_cell, result_MPPCA_multi, prms_multi, ensemble_set_A, ensemble_set_B, ensemble_set_C, ensemble_set_D)
%% comment
% 211118: 1st ver
% 220629: updated for multiple calculation
%% system parameter
shift = 1; % system parameter, don't change.
%% decompose data and parameter
result_MPPCA1 = result_MPPCA_multi.result_MPPCA1;
result_MPPCA2 = result_MPPCA_multi.result_MPPCA2;
result_MPPCA3 = result_MPPCA_multi.result_MPPCA3;
result_MPPCA4 = result_MPPCA_multi.result_MPPCA4;

binarize_thr          = prms_multi.binarize_thr; % Calcium trace binalizing threshold.
shuffle_mode          = prms_multi.shuffle_mode; % 1:randperm, 2:circshift

thr_percentile        = prms_multi.thr_percentile; % def: 1%, top-bottom %-thresholding
target_session_num    = prms_multi.target_session_num; 
ensemble_overlap_mode = prms_multi.ensemble_overlap_mode;
iteration             = prms_multi.iteration; 
%% cell sorting
temp_ca_data = result_data_cell{target_session_num+shift,7};
ensemble_cellIDs{1} = result_MPPCA1.neuron_sig_IDs{1, ensemble_set_A};
ensemble_cellIDs{2} = result_MPPCA2.neuron_sig_IDs{1, ensemble_set_B};
ensemble_cellIDs{3} = result_MPPCA3.neuron_sig_IDs{1, ensemble_set_C};
ensemble_cellIDs{4} = result_MPPCA4.neuron_sig_IDs{1, ensemble_set_D};

% ### debug, chase all cell ids
% ensemble_cellIDs_all = [ensemble_cellIDs{1}; ensemble_cellIDs{2}; ensemble_cellIDs{3}; ensemble_cellIDs{4}]; 
%% ### figure Venn diagram
% figure
% setListData = {ensemble_cellIDs{1}, ensemble_cellIDs{2}, ensemble_cellIDs{3}, ensemble_cellIDs{4}};
% setLabels = ["P1", "P2", "P3", "P4"];
% h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);
% xticklabels([]); yticklabels([]); title('Venn diagram')
%% find unique cell ids
P1_unique = setdiff(ensemble_cellIDs{1}, [ensemble_cellIDs{2}; ensemble_cellIDs{3}; ensemble_cellIDs{4}]);
P2_unique = setdiff(ensemble_cellIDs{2}, [ensemble_cellIDs{1}; ensemble_cellIDs{3}; ensemble_cellIDs{4}]);
P3_unique = setdiff(ensemble_cellIDs{3}, [ensemble_cellIDs{1}; ensemble_cellIDs{2}; ensemble_cellIDs{4}]);
P4_unique = setdiff(ensemble_cellIDs{4}, [ensemble_cellIDs{1}; ensemble_cellIDs{2}; ensemble_cellIDs{3}]);
%% Calculate connectivity
if ensemble_overlap_mode == 0
[Cells_sum] = fxn_MPPCA_cal_connect_4patterns(temp_ca_data, P1_unique, P2_unique, P3_unique, P4_unique, binarize_thr); % exclude overlap
elseif ensemble_overlap_mode == 1
[Cells_sum] = fxn_MPPCA_cal_connect_4patterns(temp_ca_data, ensemble_cellIDs{1}, ensemble_cellIDs{2}, ensemble_cellIDs{3}, ensemble_cellIDs{4}, binarize_thr); % exclude overlap
end
%% Calculate shuffled connetctivity
% shuffle data generation
% shuffle_mode = 2;
% shuffled_data = fxn_circ_shuffling(data_input,iteration);
% temp_shuffled_data = fxn_data_shuffling(temp_ca_data, iteration, shuffle_mode); % each cell shuffling
temp_ca_shuffled{1,1} = fxn_ensemble_shuffling(temp_ca_data, iteration, shuffle_mode); % ensemble shuffling
temp_ca_shuffled{2,1} = fxn_ensemble_shuffling(temp_ca_data, iteration, shuffle_mode); % ensemble shuffling
temp_ca_shuffled{3,1} = fxn_ensemble_shuffling(temp_ca_data, iteration, shuffle_mode); % ensemble shuffling
temp_ca_shuffled{4,1} = fxn_ensemble_shuffling(temp_ca_data, iteration, shuffle_mode); % ensemble shuffling
%% shuffle data cal, This section takes time. 
Cells_sum_shuffle = cell(iteration,1); % for speed up cal.
for i_shuffle = 1:iteration
[Cells_sum_shuffle{i_shuffle,1}] = fxn_MPPCA_cal_connect_4patterns_shuffle ...
    (temp_ca_shuffled, P1_unique, P2_unique, P3_unique, P4_unique, binarize_thr, i_shuffle);
end
% disp(['Gerenating surrogate data with ', num2str(iteration), ' iterations.']);
%% Distribution thresholding
% iteration = 1000;
% thr_percentile = 1; % 1%-thresholding

[Cells_sum_shuffle_res] = fxn_MPPCA_cal_distribution_thr_realine(Cells_sum_shuffle, iteration, thr_percentile);
%%
% figure('Position', [50 50 600 200]) ;
% subplot(121);
% histogram(cell2mat(Cells_AB_shuffle_index))
% hold on
% xline(Cells_AB_sum_index,'r','Actual value','LineWidth',2)
%% output results
result_multi.Cells_sum             = Cells_sum;
result_multi.Cells_sum_shuffle_res = Cells_sum_shuffle_res;
result_multi.P1_unique = P1_unique;
result_multi.P2_unique = P2_unique;
result_multi.P3_unique = P3_unique;
result_multi.P4_unique = P4_unique;
%%
end