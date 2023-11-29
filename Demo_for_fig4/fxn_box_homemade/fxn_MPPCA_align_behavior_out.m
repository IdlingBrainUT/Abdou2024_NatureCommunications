%% Alignment for MPPCA and behavioral alignment
function fxn_MPPCA_align_behavior_out(reference_session_num_to_show, session_num_to_show, result_MPPCA, results_XY, results_task, dirname1, dirname2,dirname3)
%%
ith_decompose           = result_MPPCA.ith_decompose;
r_kt                    = result_MPPCA.r_strength_ref;
r_kt_z                  = result_MPPCA.r_strength_refz;
r_kt_full               = result_MPPCA.r_strength_target;
r_kt_full_z             = result_MPPCA.r_strength_targetz;
V_sqr_scale             = result_MPPCA.neuron_weight;
V_sqr_scale_z_cutoff    = result_MPPCA.neuron_weight_sigID;
V3_sort_id_data         = result_MPPCA.neuron_sig_IDs;

r                       = result_MPPCA.r;
prms_reactivation_SD_thr= result_MPPCA.prms_reactivation_SD_thr;
result_data_cell        = result_MPPCA.result_data_cell;
data_bin_z_full         = result_MPPCA.data_bin_z_full;
data_bin_z              = result_MPPCA.data_bin_z;
V_sqr_scale_z           = result_MPPCA.V_sqr_scale_z;
negative2               = result_MPPCA.negative2;
bin_frame_num           = result_MPPCA.bin_frame_num;

%% Post processing for SeqNMF occurrence 
% cell_file1 = {};
% dir_name1 = (iSeq_out_filename);
% file1_name = (file1_name);
% dir_name2 = ('_H');
% dir_name4 = ('.csv');
% dir_name = [dir_name1, file1_name, dir_name2, dir_name4];
% csv_numeric = csvread(dir_name);
% csv_numeric = zscore(csv_numeric);
% disp('Finishied occurrence csv reading!')
%% Add time info to CS
shift = 1;
down_sample_factor = bin_frame_num; % SeqNMFのダウンサンプリングファクターと合わせる。
nVista_fs = 20;
session_num = session_num_to_show +shift;
reference_session_num = reference_session_num_to_show +shift;
    % recording_dur = size(csv_numeric,1) * down_sample_factor/nVista_fs;
MPPCA_t = result_data_cell{session_num,5};
    % csv_t = [1:size(csv_numeric,1)]/down_sample_factor;
    % csv_num = [1:size(csv_numeric,1)];
    % csv_heuristic = csv_numeric(:,pattern_ID_target);
    % csv_heuristic_null = [csv_heuristic, zeros(size((csv_heuristic),1),4)]; % this is for just merged plot. no meaning.
    % csv_schema = csv_numeric(:,pattern_ID_all)';
    % csv_schema_null = [csv_heuristic, zeros(size((csv_heuristic),1),4)]; % this is for just merged plot. no meaning.
xlim_end = MPPCA_t(end);
%% give time info

if session_num == 2
result_data_cell_t = [1 : result_data_cell{session_num,8}];
elseif session_num > 2
result_data_cell_t = [result_data_cell{session_num-1,8}+1 : result_data_cell{session_num,8}];
end
%% Post processing for spectrogram,  Concatenate directory names 
tank_dir = [dirname1, dirname2, dirname3];
fig_title = string(dirname3);
xlim_input = [0 xlim_end]; % input sec

%% Output results file
total_Whel_time_val = results_task.total_Whel_time_val;
total_Lick_time_val = results_task.total_Lick_time_val;
total_SOL_time_val  = results_task.total_SOL_time_val;
Whel_table          = results_task.Whel_table;
total_Stage_time_val= results_task.total_Stage_time_val;

if isfield(results_XY, 'XY_posi_filled_t')
    XY_t      = results_XY.XY_posi_filled_t;
    XY_x_norm = results_XY.x_norm;
    XY_y_norm = results_XY.y_norm;    
else
    XY_t      = [1:size(results_XY,1)]/20;
    XY_x_norm = results_XY(:,1);
    XY_y_norm = results_XY(:,2);   
end            
%% figure for reference_z
% reference data fig, ensemble sorting for check
% for i = 1:1 % for debug
for i = 1:r

fig_show(i) = figure('Position',[100,0,700,700]); %[left bottom width height] 


h1 = subaxis(8,1,1, 'SpacingVert',0.03, 'MR',0.05, 'ML',0.1); 
    plot(XY_t ,XY_x_norm, 'k', 'LineWidth',0.5); hold on
    plot(XY_t ,XY_y_norm, 'r', 'LineWidth',1); xlim(xlim_input); 
    ylabel(sprintf('XY \n position')); xticklabels([]);% legend('X-position','Y-position'); 
title(['Backprojected reactivation strength for MPPCA-ICA ensemble #', ...
        num2str(i), ' of ' ,num2str(r),' in ', cell2mat(result_data_cell(session_num,1))]);
    grid on; box off;
    
h2 = subaxis(8,1,2); plot(total_Whel_time_val(:,1),total_Whel_time_val(:,2),'b','LineWidth',1);    xlim(xlim_input);  % xlim second
ylabel(sprintf('Cue \n count'));xticklabels([]);grid on; box off;
h3 = subaxis(8,1,3); plot(total_Lick_time_val(:,1),total_Lick_time_val(:,2), 'k', 'LineWidth',1); xlim(xlim_input); 
ylabel(sprintf('Total \n Lick'));xticklabels([]); grid on;  box off; %ylabel('Total Lick');
h4 = subaxis(8,1,4); plot(total_SOL_time_val(:,1),total_SOL_time_val(:,2),'r', 'LineWidth',1);    xlim(xlim_input);  % xlim second
ylabel(sprintf('Correct \n Lick')); grid on; box off; % xticklabels([]);

%% セッション抜いてzする　全体zではない。
r_kt_full_z = zscore(r_kt_full(result_data_cell_t,i));

r_kt_full_z_above_SD = find(r_kt_full_z >= prms_reactivation_SD_thr);
r_kt_full_z_below_SD = find(r_kt_full_z < prms_reactivation_SD_thr);

r_kt_full_z_color_above_thr = r_kt_full_z;  
r_kt_full_z_color_above_thr(r_kt_full_z_below_SD) = NaN;
r_kt_full_z_color_above_thr_dash = ones(size(result_data_cell_t,2),1)*prms_reactivation_SD_thr;
%%
% h2 = subaxis(6,1,2,'SpacingVert',0.01) ;
h5 = subaxis(8,1,5);

plot(MPPCA_t, r_kt_full(result_data_cell_t,i),'b'); grid on;
ylabel({'Reactivation';'strength (AU)'}); xlim([0 MPPCA_t(end)]); %ylim([0 1])
xticklabels([]); hold off;  box off;

%%
% h3 = subaxis(6,1,3); 
h6 = subaxis(8,1,6:8);

plot(MPPCA_t, r_kt_full_z,'k'); hold on; grid on;
area(MPPCA_t, r_kt_full_z_color_above_thr,'FaceColor','r','EdgeColor','r'); hold on
plot(MPPCA_t, r_kt_full_z_color_above_thr_dash, 'r--'); 
ylabel({'Reactivation';'strength (SD)'}); xlim([0 MPPCA_t(end)]); 
 xlabel('Time (s)'); hold off;  box off;
%% linkaxes
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 9, 'FontName','Arial');
linkaxes([h1,h2,h3,h4,h5,h6], 'x');
%%
end
%% plot XY
figure('Position', [500 100 150 450]);

subplot(311); plot(XY_x_norm, XY_y_norm, 'ko', 'MarkerSize',0.5);
xlabel('Position X'); ylabel('Position Y')    

subplot(312); plot(total_SOL_time_val(:,1),total_SOL_time_val(:,3) );
xlabel('# of Drink'); ylabel('Time(s)')

subplot(313); 
histogram(Whel_table(:,3), 60, 'BinLimits',[0,60],'Normalization','probability');
xlabel('Run time (s)'); ylabel('Frequency')

% subplot(414); plot(XY_t,XY_x_norm, 'k', 'LineWidth',0.5); 
%     hold on
%               plot(XY_t,XY_y_norm, 'r', 'LineWidth',1); 
% legend('X-position','Y-position'); xlim(xlim_input);
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 10, 'FontName','Arial');
%% save fig
fig_show_name = ['result_MPPCA_aligned_fig', ...
    '_ref', num2str(reference_session_num_to_show), '_', cell2mat(result_data_cell(reference_session_num,1)), ...
    '_tar', num2str(session_num_to_show),   '_', cell2mat(result_data_cell(session_num,1)),'_', num2str(r),'patterns_reference'];
savefig(fig_show, fig_show_name, 'compact')
%%
end