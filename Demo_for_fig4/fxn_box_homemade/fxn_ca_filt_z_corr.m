function [ca_filt_z] = fxn_ca_filt_z_corr_figure(results, total_session_num)
%%
% C = cat(dim,A,B)
cell_temp = results; % for debug
ca_filt_z = [];
for i = 1:total_session_num
temp_z=cell_temp{i,3}; ca_filt_z=cat(1,ca_filt_z,temp_z); 
end

% ca_filt_z = ca_filt_data;

corr_t_res = corrcoef(ca_filt_z'); corr_t_scale = [1:size(corr_t_res,1)];
corr_n_res = corrcoef(ca_filt_z); corr_n_scale = [1:size(corr_n_res,1)];
%%
figure('Position',[20 20 1000 800]);
% subplot(211); imagesc(ca_filt_z'); clim([0 3]);subplot(212); 
imagesc(corr_t_scale, corr_t_scale, corr_t_res); clim([0 0.2]);colormap('parula'); colorbar; clim([-0.3 0.3]);
xticks([0:200:10000]); yticks([0:200:10000]); xlabel('Time(s)');ylabel('Time(s)'); grid on
ax = gca; ax.TickDir = 'both';

% figure('Position',[40 40 700 600]);
% imagesc(corr_n_scale, corr_n_scale, corr_n_res); ;colormap('parula'); colorbar; clim([0 0.2]);
% xticks([0:10:700]); yticks([0:10:700]); xlabel('Neuron#');ylabel('Neuron#'); grid on
% ax = gca; ax.TickDir = 'both';

end
%%