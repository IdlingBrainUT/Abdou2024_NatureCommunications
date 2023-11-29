function [ca_filt_data] = fxn_butter_filt_cutoff(ca_raw_data, sample_fps, highpass, factor, cutoff_val)
%% for debug,
% sample_fps   = 20; % 20hz
% highpass     = 0.01; % 30 sec highpass filter
% factor       = 1; % factor
%% filtering code
[b,a] = butter(factor, highpass /(sample_fps/2),'high'); %
% [b,a] = cheby2(factor, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
ca_filt_data = filtfilt(b, a, ca_raw_data); %ゼロ位相filtfiltで ハイパス　フィルターを通す  
%% cutoff section
reply = input('Would you like to do cutoff with set value? (y/n): ','s');
% reply = 'y'; % skipping

if strcmp(reply,'y') % 
disp1 = ('You chose YES. cutoff done!');
disp([disp1, ' with cutoff value ', num2str(cutoff_val), '.']) 
negative = find(ca_filt_data < cutoff_val); % 2SD以下をカット
ca_filt_data(negative) = zeros(size(negative));
display('Finish butter filt and cutoff!')
    
elseif  strcmp(reply,'n') % 
disp('You chose NO. Only butter filt done!')
 
else
error('Unexpected situation')
end
%% figure
figure;
subplot(311); imagesc(ca_raw_data'); caxis([0 3]); title('Before')
subplot(312); imagesc(ca_filt_data'); caxis([0 3]); title('After')
subplot(313); plot(mean(ca_raw_data', 1),'b'); hold on;
              plot(mean(ca_filt_data', 1), 'r');
%%    
end