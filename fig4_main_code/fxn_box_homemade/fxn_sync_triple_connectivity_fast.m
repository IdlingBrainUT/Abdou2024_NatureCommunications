function [raster_triple_connectivity_t_index, raster_triple_connectivity_index] = ...
    fxn_sync_triple_connectivity_fast(data_input_temp, Cell_A_ID, Cell_B_ID, Cell_C_ID)
%% update info
% 211111: Calculation section for combination matrix is improved. 
%% for debug
% ### mymap_gradation_set, 1->Black, 2-> blue and green, 3-> red and green
% data_input_temp = G1_results{1,6};
% Cell_type_A_ID = G1_results{2,5};
% 
% Cell_A_ID = Cell_type_ABC_pure_A_ID;
% Cell_B_ID = Cell_type_ABC_pure_B_ID;
% Cell_C_ID = Cell_type_ABC_pure_C_ID;
% mymap_gradation_set = 2;
%% data load
cal_range_for_bin_data = 1:(size(data_input_temp,1));
%% Occurance calculation
%     data_input_temp_sum = sum(data_input_temp,2);
%     [hist_B1, ~] = find(data_input_temp_sum >0);
%     data_for_hist1 = data_input_temp_sum(hist_B1);
%     
%     figure; subplot(223);
%      hist = histogram(data_for_hist1(:,1));
%    
% %     hist.NumBins = 100;
%     hist.BinEdges = [0:2:100];
% %     ylim([0 20])
%     
%     title_1 = ('cut-frame-num=');
%     title_2 = (', cal-range=');
%     title_3 = (' in sec');
% %     title_123 = [title_1, num2str(cut_frame_num), title_2, num2str(xlim_range), title_3];
%     
% %     title(title_123)
%     
%          xlabel('# of Synchronized cells','FontSize',12,'FontWeight','bold','Color','k'); 
%          ylabel('Occurance','FontSize',12,'FontWeight','bold','Color','k'); 
% 
% %          set(gca,'xtick', 1:2:100)
% %          set(gca,'xticklabel',(1:50))
%          xtickangle(90);
%% figure multi color raster
% % mymap_gradation = [1 1 1;0 0 1; 1 0 0; 0 0.5 0]; % white111, blue001, red100, deepgreen0 0.5 0
% % ax4 = subplot(224); imagesc(raster_connection_sum'); colormap(ax4,mymap); grid on; 
% % figure;
% clim_range = [0 1];
% xlim_range = [0 size(data_input_temp,1)];
% bin_frame_num = 1; % already binned
% 
% % ID check
% posi_ID = [Cell_A_ID, Cell_B_ID, Cell_C_ID];
% nega_ID = (1:size(data_input_temp,2));
% nega_ID((posi_ID)) = []; 
% 
% data_color1 = data_input_temp(:,nega_ID);
% data_color2 = (data_input_temp(:,Cell_A_ID)*2);
% data_color3 = (data_input_temp(:,Cell_B_ID)*3);
% data_color4 = (data_input_temp(:,Cell_C_ID)*4);
% 
% data_color_all = [data_color1,data_color2,data_color3,data_color4];
% 
% % figure;
% mymap_4colors = [1 1 1; 0 0 0; 0 0 1; 1 0 0; 0 0.5 0]; 
% % (0)white111, (1)black000; (2)blue001, (3)red100, (4)deepgreen0 0.5 0
% ax1 = subplot(221); 
% imagesc((1:size(data_color_all,1)), (1:size(data_color_all,2)) ,data_color_all')
% xlim(xlim_range); %clim(clim_range);
% %          title(title_123)
% 
%          % bold off
%          xlabel('Time (s)','FontSize',12,'Color','k'); 
%          ylabel('Neuron#','FontSize',12,'Color','k'); 
% 
% %          xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
% %          ylabel('Neuron#','FontSize',12,'FontWeight','bold','Color','k'); 
% %          
%          colormap(ax1,mymap_4colors); grid on; %xticks([xticks_range]); 
%% for debug, Input cell ID
% Cell_type_A_ID = [1,2,3,4,5];
% Cell_type_B_ID = [6,7,8,9,10,11,12,13,14,15];
% Cell_type_C_ID = [31,32,33];
%% 
Cells_A = data_input_temp(:,Cell_A_ID);
Cells_B = data_input_temp(:,Cell_B_ID);
Cells_C = data_input_temp(:,Cell_C_ID);
%% 211111 added. Sum-methond

Cells_A_sum = sum(Cells_A, 2);
Cells_B_sum = sum(Cells_B, 2);
Cells_C_sum = sum(Cells_C, 2);

for i = 1:size(Cells_A_sum,1)
Cells_ABC_sum_res(i,1) = Cells_A_sum(i)*Cells_B_sum(i)*Cells_C_sum(i);
end
%% raster_normalization new
raster_triple_connectivity_t_index = Cells_ABC_sum_res/ (size(Cells_A,2)*size(Cells_B,2)*size(Cells_C,2));
raster_triple_connectivity_index = sum(raster_triple_connectivity_t_index)/size(Cells_A,1);

%%
% %% Calculation triple interaction
% for i = 1:size(Cells_A,2)
%     for ii=1:size(Cells_B,2) 
%        for iii=1:size(Cells_C,2)
%             for iiii = 1:size(Cells_A,1)
% 
%                     if Cells_A(iiii,i)==1 && Cells_B(iiii,ii)==1 && Cells_C(iiii,iii)==1
%                         Cells_A_B_C(i,ii,iii,iiii) =1; 
%                     else
%                         Cells_A_B_C(i,ii,iii,iiii) =0;
%                     end
%             end
%         end
%     end
% end
% %% fix code
% if exist('Cells_A_B_C')
%     disp('Cells_A_B_C is normally calculated.')   
% else
%     disp('Set Cells_A_B_C as 4dim-zeros matrix.')   
%     Cells_A_B_C = zeros(size(Cells_A,2),size(Cells_B,2),size(Cells_C,2),size(Cells_A,1));
%     i = 0;
%     ii = 0;
%     iii = 0;
% end
% %% Dimension squeeze
% temp_A = sum(Cells_A_B_C,3); 
% temp_A_sqz = squeeze(temp_A);
% temp_A2 = sum(temp_A_sqz,2); 
% temp_A_sqz2 = squeeze(temp_A2);
% temp_A_sqz_logical = temp_A_sqz2'> 0;
% temp_A_sqz_logical_sum = sum(temp_A_sqz_logical,2);
% 
% temp_B = sum(Cells_A_B_C,3); 
% temp_B_sqz = squeeze(temp_B);
% temp_B2 = sum(temp_B_sqz,1); 
% temp_B_sqz2 = squeeze(temp_B2);
% temp_B_sqz_logical = temp_B_sqz2'> 0;
% 
% temp_C = sum(Cells_A_B_C,1); 
% temp_C_sqz = squeeze(temp_C);
% temp_C2 = sum(temp_C_sqz,1); 
% temp_C_sqz2 = squeeze(temp_C2);
% temp_C_sqz_logical = temp_C_sqz2'> 0;
% %% Calculate entire occurence
% Cells_ABC_res_occurence_raster_fig = [temp_A_sqz_logical*1, temp_B_sqz_logical*2, temp_C_sqz_logical*3];
% Cells_ABC_res_occurence_sum = sum(Cells_ABC_res_occurence_raster_fig,2);
% % Cells_A_B_C_res_occurence_vec = Cells_A_B_C_res_occurence_sum > 0;
% %% ### triple calculation end
% 
% %% fix code
% if exist('Cells_A_B_C')
%     disp('Cells_A_B_C is normally calculated.')   
% %% Figure synchronized raster
% 
%     title_0 = ('Sorted ');
% %     title_0123 = [title_0, title_1, num2str(cut_frame_num), title_2, num2str(xlim_range), title_3];
% 
% % figure; 
% % ax4 = subplot(224); imagesc(raster_connection_sum'); colormap(ax4,mymap); grid on; 
% ax2 = subplot(222); 
% imagesc( (1:size(Cells_ABC_res_occurence_raster_fig',2))*bin_frame_num , (1:size(Cells_ABC_res_occurence_raster_fig',1)) ,Cells_ABC_res_occurence_raster_fig')
% % xlim(xlim_range); clim(clim_range);
% %          title(title_0123)
% 
%          % bold off
%          xlabel('Time (s)','FontSize',12,'Color','k'); 
%          ylabel('Neuron#','FontSize',12,'Color','k'); 
%          
% %          xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
% %          ylabel('Neuron#','FontSize',12,'FontWeight','bold','Color','k'); 
% %          
% %         colorbar
% mymap_gradation_set = 2;
% 
%         if mymap_gradation_set == 0;  mymap_gradation = [1 1 1; 0 0 0];
%     elseif mymap_gradation_set == 1;  mymap_gradation = [1 1 1; 0 0 1; 0 0.5 0];
%     elseif mymap_gradation_set == 2;  mymap_gradation = [1 1 1; 0 0 1; 1 0 0; 0 0.5 0];
%         else display('Input correct value!')
%         end
%    
% % mymap_gradation = [1 1 1;0 0 1; 1 0 0; 0 0.5 0]; % white111, blue001, red100, deepgreen0 0.5 0
% colormap(ax2,mymap_gradation); grid on; %xticks([xticks_range]);   
% 
% % 色	RGB 3 成分
% % 黄	[1 1 0]
% % マゼンタ	[1 0 1]
% % シアン	[0 1 1]
% % 赤	[1 0 0]
% % 緑	[0 1 0]
% % 青	[0 0 1]
% % 白	[1 1 1]
% % 黒	[0 0 0]
% %% fixing code
% % if exist('Cells_A')
% %     disp('Cells_A is normally calculated.')   
% % else
% %      Cells_A = zeros(10,2); % 適当に10を代入 10 x 10 x 10 matrix
% %     disp('Cells_A is set zero.')
% % end
% % 
% % if exist('Cells_B')
% %     disp('Cells_B is normally calculated.')   
% % else
% %      Cells_B = zeros(10,2); % 適当に10を代入 10 x 10 x 10 matrix
% %     disp('Cells_B is set zero.')
% % end
% %% Calculate synchronized connectivity 
% 
% raster_A = temp_A_sqz_logical;
% raster_B = temp_B_sqz_logical;
% raster_C = temp_C_sqz_logical;
% 
% % i_A  = 1:size(raster_A,1);
% % ii_A = 1:size(raster_A,2);
% 
% % i_B = 1:size(raster_B,1);
% % ii_B = 1:size(raster_B,2);
% 
% raster_connection = zeros(size(raster_A,1),size(raster_A,2),size(raster_B,2)); % initialize
% 
% % for i_A = 1:size(raster_A,1)    
% %     for ii_A  = 1:size(raster_A,2)    
% %         for ii_B = 1:size(raster_B,2)
% %             
% %            if   Cells_A(i_A,ii_A)==1 && Cells_B(i_A, ii_B) ==1
% %                raster_connection(i_A,ii_A,ii_B) = 1;
% %            else
% %                raster_connection(i_A,ii_A,ii_B) = 0;
% %            end
% %            
% %         end
% %     end
% % end
% 
% for  iiii = 1:size(raster_A,1)
%     for i  = 1:size(raster_A,2)    
%         for ii = 1:size(raster_B,2)
%             for iii = 1:size(raster_C,2)
%                 
%            if   raster_A(iiii,i)==1 && raster_B(iiii,ii) ==1 && raster_C(iiii,iii) ==1
%                raster_connection(iiii,i,ii,iii) = 1;
%            else
%                raster_connection(iiii,i,ii,iii) = 0;
%            end
%             
%             end
%         end
%     end
% end
% 
% raster_connection_sum_3d = sum(raster_connection,4);
% raster_connection_sum_2d = sum(raster_connection_sum_3d,3);
% raster_connection_sum_1d = sum(raster_connection_sum_2d,2);
% %% raster_normalization old
% raster_triple_connectivity_t_index = raster_connection_sum_1d/ (size(raster_A,2)*size(raster_B,2)*size(raster_C,2));
% raster_triple_connectivity_index = sum(raster_triple_connectivity_t_index)/size(raster_A,1);
% 
% else
%     i = 0;
%     ii = 0;
%     iii = 0;
%     raster_triple_connectivity_t_index = zeros(size(data_input_temp,1),1);
%     raster_triple_connectivity_index = 0; 
%     disp('Cells_ABC is set zero.')
% end
% 
% raster_connectivity_index_str = num2str(raster_triple_connectivity_index);
% disp_index = ['Occurence value is ', raster_connectivity_index_str];
% disp(disp_index);
%     
% ax4 = subplot(224); plot(raster_triple_connectivity_t_index'); grid on; %colormap(ax4,mymap); 
% % ax1 = subplot(1,2,1); colormap(ax1,parula)
% xlabel('Time (s)','FontSize',12,'FontWeight','bold','Color','k'); 
% ylabel('Connectivity','FontSize',12,'FontWeight','bold','Color','k'); 
% %% あとで追加
%    title_1 = ('cut-frame-num=');
%     title_2 = (', cal-range=');
%     title_3 = (' in sec');
% %     title_123 = [title_1, num2str(cut_frame_num), title_2, num2str(xlim_range), title_3];
% %     title(title_123)
% 
% %% raster_connection check
% % figure; 
% % subplot(311); imagesc(raster_connection_sum'); 
% % subplot(312); imagesc(raster_connection_sum_1d')
% % subplot(313); imagesc(raster_connectivity_t_index')
% 
% %% Occurance calculation for histogram
% 
% %     data_range_sum = Cells_A_B_res_occurance_sum(cal_range_for_bin_data,:);
% %     [hist_B2, ~] = find(data_range_sum>0);
% %     data_for_hist2 = data_range_sum(hist_B2);
% % 
% % %     figure; 
% %    subplot(223)
% %    hist_sort = histogram(data_for_hist2);%
% % %     hist.NumBins = 100;
% % %     hist.BinEdges = [0:2:100];
% % %     ylim([0 20])
% %     
% %     title_1 = ('cut-frame-num=');
% %     title_2 = (', cal-range=');
% %     title_3 = (' in sec');
% % %     title_123 = [title_1, num2str(cut_frame_num), title_2, num2str(xlim_range), title_3];
% %     
% % %     title(title_123)
% %     
% %          % bold off
% %          xlabel('# of Synchronized cells','FontSize',12,'Color','k'); 
% %          ylabel('Occurance','FontSize',12,'Color','k'); 
% % 
% % 
% % %          xlabel('# of Synchronized cells','FontSize',12,'FontWeight','bold','Color','k'); 
% % %          ylabel('Occurance','FontSize',12,'FontWeight','bold','Color','k'); 
% % 
% % %          set(gca,'xtick', 1:2:100)
% % %          set(gca,'xticklabel',(1:50))
% %          xtickangle(90);
% %          
%          set(findobj(gcf, 'Type', 'Axes'), 'FontSize', 8, 'FontName','Arial'); %grid on;
% % 
% % %          disp('show histogram results in sorted data'); 
% %%
% 
end
