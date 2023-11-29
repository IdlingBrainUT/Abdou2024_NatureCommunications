%% Calculate correlation matrix
function  [res_corr_mat_t, res_corr_mat_n, data_cell_mean] = fxn_ca_correlation_matrix_posi_only(data_cell, total_session_num, Cell_IDs_posi, ref_sec_num_from_end)
%% for debug
cell_temp = data_cell;
shift =1; 
% ref_sec_num_from_end = 1; % for reference section number frmo end
%%
% C = cat(dim,A,B)
res_time_mat = [];
for i = 1:total_session_num
temp_z = cell_temp{i+shift,7}(:,[Cell_IDs_posi]); % skip ;Cell_IDs_nega
res_time_mat = cat(1,res_time_mat,temp_z); 
end

res_corr_mat_t = corrcoef(res_time_mat'); corr_t_scale = [1:size(res_corr_mat_t,1)];
res_corr_mat_n = corrcoef(res_time_mat); corr_n_scale = [1:size(res_corr_mat_n,1)];
%% figure for time correlation
figure('Position',[20 20 1000 800]);
imagesc(corr_t_scale, corr_t_scale, res_corr_mat_t);
xticks([1; cell2mat(data_cell(1+shift:total_session_num+1,8))]); 
xticklabels(data_cell(1+shift:total_session_num+1,1)); 
yticks([1; cell2mat(data_cell(1+shift:total_session_num+1,8))]); 
yticklabels(data_cell(1+shift:total_session_num+1,1)); 

colormap(fxn_redblue); colorbar; clim([-0.2 0.2]);
xlabel('Time(s)');ylabel('Time(s)'); grid on
c = colorbar; c.Label.String = 'Correlation index';
ax = gca; ax.TickDir = 'both'; ax.XAxisLocation = 'top';
set(gca, 'FontSize',12, 'FontName','Arial'); set(c, 'FontSize',12, 'FontName','Arial'); 
%% figure for cell-cell correlation
figure('Position',[40 40 1000 800]);
imagesc(corr_n_scale, corr_n_scale, res_corr_mat_n); colormap(parula); colorbar; clim([-0.3 0.3]);
xticks([0:20:700]);
yticks([0:20:700]); 

colormap(fxn_redblue); colorbar; clim([-0.2 0.2]);
xlabel('Neuron#');ylabel('Neuron#'); grid on
c = colorbar; c.Label.String = 'Correlation index'; 
ax = gca; ax.TickDir = 'both'; ax.XAxisLocation = 'top';
set(gca, 'FontSize',12, 'FontName','Arial'); 
%% Calculate mean correlation value
for i = 1:total_session_num
    if i == 1
data_cell{i+shift,9}  = (1: data_cell{i+shift,8});
data_cell{i+shift,10} = res_corr_mat_t(:, data_cell{i+shift,9});  %  cumulative_time; vectorize   
    end
    if i > 1
data_cell{i+shift,9}  = ((1+data_cell{i,8}): data_cell{i+shift,8}); 
data_cell{i+shift,10} = res_corr_mat_t(:, (1+data_cell{i,8}): data_cell{i+shift,8});  %  cumulative_time; vectorize
    end
data_cell{i+shift,11} = mean(data_cell{i+shift,10},2); 
end
for ii = 1:total_session_num
 for iii =  1:total_session_num
data_cell_mean(iii,ii)  = mean(data_cell{ii+shift,11}(data_cell{iii+shift,9}),1);
 end
end
%
%% figure simpel mean correlation matrix
figure('Position', [100, 100, 800, 300]);  
subaxis(1,2,1, 'SpacingVert',0.02, 'MR',0.05, 'ML',0.1, 'Padding', 0.05); 
imagesc(data_cell_mean);
xticklabels(data_cell(1+shift:total_session_num+1,1)); 
yticklabels(data_cell(1+shift:total_session_num+1,1));
xticks([1:total_session_num+1]); yticks([1:total_session_num+1]); xtickangle(45); title('Correlatin matrix')
set(gca, 'FontSize',10, 'FontName','Arial'); %colormap(fxn_redblue); colorbar; clim([-0.015 0.015]);

subaxis(1,2,2)
plot(data_cell_mean(:,end- ref_sec_num_from_end), 'g', 'LineWidth',2)
xticklabels(data_cell(1+shift:total_session_num+1,1)); 
xticks([1:total_session_num+1]); xtickangle(45); ylabel('Correltation mean'); title('Correlation (vs. refernce)')
set(gca, 'FontSize',10, 'FontName','Arial'); 
%%
end
%%