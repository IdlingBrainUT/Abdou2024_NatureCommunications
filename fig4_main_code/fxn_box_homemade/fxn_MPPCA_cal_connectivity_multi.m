function [result_final_cell, result_multi] = fxn_USM_MPPCA_cal_connectivity_multi(result_data_cell, result_MPPCA_multi, prms_multi)
%% comment
% 211118: 1st version
%% (Multi mode) Parameter input section
% ### input saved filename. To protect calculation resutls.
dt = datetime('now');
filename_DateString = datestr(dt,'yymmdd_HHMM');
filename_for_save = (['Res_multi_cal_',filename_DateString,'.mat']); 
%% decompose prms
target_session_num = prms_multi.target_session_num ;
binarize_thr       = prms_multi.binarize_thr       ;
shuffle_mode       = prms_multi.shuffle_mode       ;
thr_percentile     = prms_multi.thr_percentile     ;
iteration          = prms_multi.iteration          ;                    
ensemble_set_A     = prms_multi.ensemble_set_A     ;
ensemble_set_B     = prms_multi.ensemble_set_B     ;
ensemble_set_C     = prms_multi.ensemble_set_C     ;
ensemble_set_D     = prms_multi.ensemble_set_D     ;
%%
%% CPU cal. ver for multi-wise connectivity
% fxn_MPPCA_cal_connectivity_multi_core(result_data_cell, result_MPPCA_multi, prms_multi, ensemble_set_A, ensemble_set_B, ensemble_set_C, ensemble_set_D)
tic;
result_multi = cell(ensemble_set_A, ensemble_set_B, ensemble_set_C, ensemble_set_D);
loop_num_current = 0; % initialize

parfor i = 1:ensemble_set_A % parfor ok
% for i = 1:ensemble_set_A % debug
    for ii = 1:ensemble_set_B
        for iii = 1:ensemble_set_C
            for iiii = 1:ensemble_set_D

% normal cpu cal.
[result_multi{i,ii,iii,iiii}] = fxn_MPPCA_cal_connectivity_multi_core(result_data_cell, result_MPPCA_multi, prms_multi, i, ii, iii, iiii); % core cal

% gpu cal.
% result_data_cell = gpuArray(result_data_cell);
% result_MPPCA_multi = gpuArray(result_MPPCA_multi);
% prms_multi = gpuArray(prms_multi);
% i = gpuArray(i);
% ii = gpuArray(ii);
% iii = gpuArray(iii);
% iiii = gpuArray(iiii);
% [gpu_temp] = arrayfun(@fxn_MPPCA_cal_connectivity_multi_core, result_data_cell, result_MPPCA_multi, prms_multi, i, ii, iii, iiii);
%  result_multi{i,ii,iii,iiii} = gather(gpu_temp);

        loop_num_current = ensemble_set_B*ensemble_set_C*ensemble_set_D*(i-1) + ensemble_set_C*ensemble_set_D*(ii-1) + ensemble_set_D*(iii-1) + iiii;
        loop_num_total   = ensemble_set_A*ensemble_set_B*ensemble_set_C*ensemble_set_D;
        disp(['Cal. is now ', num2str(loop_num_current),'/',num2str(loop_num_total),', for ',num2str(i),'/',num2str(ensemble_set_A),' - ',num2str(ii),'/',num2str(ensemble_set_B),' - ', ...
                             num2str(iii),'/',num2str(ensemble_set_C),' - ',num2str(iiii),'/',num2str(ensemble_set_D), ...
                             ' with ', num2str(iteration),' iterations.'])
            end
        end
    end
end
toc;
%% debut for checking results multi
% check1 = result_multi{1,1,1,1};
% check2 = result_multi{1,1,1,2};
%% cal ABCD index
result_cell_4D = {};
result_cell_4D_surrogate = {};

for i_res = 1:11
    for i = 1:ensemble_set_A % parfor ok
        for ii = 1:ensemble_set_B
            for iii = 1:ensemble_set_C
                for iiii = 1:ensemble_set_D
                    result_cell_4D_temp = result_multi{i,ii,iii,iiii}.Cells_sum{i_res, 4};
                    result_cell_4D_temp(isnan(result_cell_4D_temp)) = 0;
                    result_cell_4D{i_res,1}(i,ii,iii,iiii) = result_cell_4D_temp;
                    
                    result_cell_4D_surrogate_temp = ... 
                       mean(result_multi{i,ii,iii,iiii}.Cells_sum_shuffle_res.Cells_sum_shuffle_index{i_res, 1});
                    result_cell_4D_surrogate_temp(isnan(result_cell_4D_surrogate_temp)) = 0;
                    result_cell_4D_surrogate{i_res,1}(i,ii,iii,iiii) = result_cell_4D_surrogate_temp;
                end
            end
        end
    end
result_final_cell{i_res,2} = mean(cell2mat(result_cell_4D(i_res)),[1 2 3 4]);
result_final_cell{i_res,3} = mean(cell2mat(result_cell_4D_surrogate(i_res)),[1 2 3 4]);
end
%% Output results
result_final_cell(1:11,1)= {'AB','AC','AD','BC','BD','CD','ABC','ABD','ACD','BCD','ABCD'};
disp('finish calculaiton for multi');

ax1 = figure('Position',[100 100 1200 300]);
subplot(131)
plot(cell2mat(result_final_cell(:,2)),'k'); hold on 
plot(cell2mat(result_final_cell(:,3)),'--r'); hold off
xticks(1:11); xticklabels(result_final_cell(1:11,1))
title('Multiple synchronized activity')
ylabel('Synchronized index')
xlabel('Combination of patterns')
legend('Real', 'Surrogate'); box off
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');

subplot(132)
plot(cell2mat(result_final_cell(1:6,2)),'k'); hold on 
plot(cell2mat(result_final_cell(1:6,3)),'--r'); hold off
xticks(1:6); xticklabels(result_final_cell(1:6,1))
title('Multiple synchronized activity')
ylabel('Synchronized index')
xlabel('Combination of patterns')
legend('Real', 'Surrogate'); box off
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');

subplot(133)
plot(cell2mat(result_final_cell(7:11,2)),'k'); hold on 
plot(cell2mat(result_final_cell(7:11,3)),'--r'); hold off 
xticks(1:5); xticklabels(result_final_cell(7:11,1))
title('Multiple synchronized activity')
ylabel('Synchronized index')
xlabel('Combination of patterns')
legend('Real', 'Surrogate'); box off
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial');
%%
% for mat file
save(filename_for_save, 'result_final_cell', 'result_multi', '-mat');

% for pdf file
extension = ('pdf');
file_name2 = [filename_for_save, 'Figure_.', extension];
exportgraphics(ax1, file_name2, 'ContentType','vector');
%%
toc;
%%
end