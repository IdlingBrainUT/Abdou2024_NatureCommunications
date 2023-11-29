%%
function fxn_figure_task_tank_occrrence(iSeq_out_filename, file1_name, pattern_ID_target, pattern_ID_all, results_XY, results_task, dirname1, dirname2,dirname3)

%% Post processing for SeqNMF occurrence 
cell_file1 = {};
dir_name1 = (iSeq_out_filename);
file1_name = (file1_name);
dir_name2 = ('_H');
dir_name4 = ('.csv');
dir_name = [dir_name1, file1_name, dir_name2, dir_name4];
csv_numeric = csvread(dir_name);
csv_numeric = zscore(csv_numeric);
disp('Finishied occurrence csv reading!')
%% Add time info to CS
down_sample_factor = 5; % SeqNMFのダウンサンプリングファクターと合わせる。
nVista_fs = 20;
recording_dur = size(csv_numeric,1)*down_sample_factor/nVista_fs;
csv_t = [1:size(csv_numeric,1)]/down_sample_factor;
csv_num = [1:size(csv_numeric,1)];
csv_heuristic = csv_numeric(:,pattern_ID_target);
csv_heuristic_null = [csv_heuristic, zeros(size((csv_heuristic),1),4)]; % this is for just merged plot. no meaning.
csv_schema = csv_numeric(:,pattern_ID_all)';
csv_schema_null = [csv_heuristic, zeros(size((csv_heuristic),1),4)]; % this is for just merged plot. no meaning.
xlim_end = csv_t(end);
%% Post processing for spectrogram,  Concatenate directory names 
tank_dir = [dirname1, dirname2, dirname3];
fig_title = string(dirname3);
xlim_input = [0 xlim_end]; % input sec

% [results_task] = fxn_TANK_task_export(tank_dir, fig_title, xlim_input);
%% Output results file
total_Whel_time_val = results_task.total_Whel_time_val;
total_Lick_time_val = results_task.total_Lick_time_val;
total_SOL_time_val  = results_task.total_SOL_time_val;
Whel_table          = results_task.Whel_table;
total_Stage_time_val= results_task.total_Stage_time_val;

XY_t      = results_XY.XY_posi_filled_t;
XY_x_norm = results_XY.x_norm;
XY_y_norm = results_XY.y_norm;
%% skipped. figure merged plot
% fig_export1 = figure('Position',[100,100,400,400]); %[left bottom width height] 
% 
% subplot(4,1,1); imagesc(csv_t, 1:size(csv_heuristic_null,2), csv_heuristic_null');  xlim(xlim_input);
% ylabel({'Sara';'occurrence';'pattern#'}); yticklabels(string(pattern_ID_target)); clim([0 3]); 
% yticks([1:size(csv_heuristic_null,2)]); 
% hold on
% subplot(4,1,1); plot(XY_t ,XY_x_norm*2+3.5, 'w', 'LineWidth',0.5); hold on
%                 plot(XY_t ,XY_y_norm*2+3.5, 'r', 'LineWidth',1); xlim(xlim_input); % legend('X-position','Y-position'); 
%                 
% subplot(4,1,2); imagesc(csv_t, 1:size(csv_heuristic_null,2), csv_heuristic_null');  xlim(xlim_input);
% ylabel({'Sara';'occurrence';'pattern#'}); yticklabels(string(pattern_ID_target)); clim([0 3]); 
% yticks([1:size(csv_heuristic_null,2)]); 
% hold on
% subplot(4,1,2); plot(XY_t ,XY_x_norm*2+3.5, 'w', 'LineWidth',0.5); hold on
%                 plot(XY_t ,XY_y_norm*2+3.5, 'r', 'LineWidth',1); xlim(xlim_input); % legend('X-position','Y-position');            
%% figure for task tank

% xlim_input = [1500 1800]; % if ON, manual input in sec ymaze-d1
% xlim_input = [0 1500]; % if ON, manual input in sec ymaze-d1

fig_export1 = figure('Position',[100,100,700,650]); %[left bottom width height] 

title(fig_title);
% subaxis(20,1,1, 'SpacingVert',0.005, 'MR',0.05); 
h1 = subaxis(8,1,1, 'SpacingVert',0.02, 'MR',0.05, 'ML',0.1); 
                plot(XY_t ,XY_x_norm, 'k', 'LineWidth',0.5); hold on
                plot(XY_t ,XY_y_norm, 'r', 'LineWidth',1); xlim(xlim_input); xticklabels([]);% legend('X-position','Y-position'); 
h2 = subaxis(8,1,2); plot(total_Whel_time_val(:,1),total_Whel_time_val(:,2),'b','LineWidth',1);    xlim(xlim_input);  % xlim second
ylabel(sprintf('Cue \n count'));xticklabels([]);
h3 = subaxis(8,1,3); plot(total_Stage_time_val(:,1),total_Stage_time_val(:,2),'b', 'LineWidth',1); xlim(xlim_input);  % xlim second
ylabel(sprintf('Reward \n stage'));xticklabels([]);% ylabel('Reward stage');
h4 = subaxis(8,1,4); plot(total_Lick_time_val(:,1),total_Lick_time_val(:,2), 'k', 'LineWidth',1); xlim(xlim_input); 
ylabel(sprintf('Total \n Lick'));xticklabels([]);%ylabel('Total Lick');
h5 = subaxis(8,1,5); plot(total_SOL_time_val(:,1),total_SOL_time_val(:,2),'r', 'LineWidth',1);    xlim(xlim_input);  % xlim second
ylabel(sprintf('Correct \n Lick')); xticklabels([]);

h6 = subaxis(8,1,6); imagesc(csv_t, [1:numel(pattern_ID_target)], csv_numeric(:,pattern_ID_target)');  xlim(xlim_input);
ylabel({'REM+';'pattern#'}); yticklabels(string(pattern_ID_target)); clim([0 3]);yticks([1:numel(pattern_ID_target)]);xticklabels([]);

h7 = subaxis(8,1,7:8); imagesc(csv_t, [1:numel(pattern_ID_all)], csv_numeric(:,pattern_ID_all)'); xlim(xlim_input); colormap jet
ylabel({'ALL';'pattern#'}); xlabel('Time (s)'); yticklabels(string(pattern_ID_all));  yticks([1:numel(pattern_ID_all)]) ;clim([0 3]);

grid on
set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 10, 'FontName','Arial');
colormap('parula');
% set(gca, 'FontSize',12, 'FontName','Arial'); colormap('parula');

% #### meta data .emp file saving
% saveas(fig_export1,'fig_spectro1','meta')
linkaxes([h1,h2,h3,h4,h5,h6,h7], 'x');

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
%% REM_frame_Occ_ranksum
% occ_temp = csv_numeric;
% REM_range = [200 230; 750 830];
% upsampling_factor = 5; % 1 frame of iseq occurence takes 200ms. 200ms x 5 = 1 second.

% occ_temp = csv_numeric;
% upsampling_factor = 5; % 1 frame of iseq occurence takes 200ms. 200ms x 5 = 1 second.
% 
% [occ_ranksum_res] = fxn_Occ_ranksum(occ_temp, REM_range_for_stat, upsampling_factor); %  REM_range = [200 230; 750 830];
%%
%%
end

