%% correlation matrix

input_data = rand(20);

data_temp = input_data;  % input data

corr_res = corrcoef(data_temp);
corr_t = [1:size(corr_res,1)];

figure('Position',[20 20 400 400]);
subplot(121)
imagesc(corr_t, corr_t, corr_res) 

subplot(122);
corr_total = sum(corr_res,1);
plot(corr_total);

%%