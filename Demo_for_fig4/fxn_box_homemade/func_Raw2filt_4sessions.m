function func_Raw2filt_4sessions(ca_raw_data, save_mat_name, s1_range, s2_range, s3_range, s4_range)

%% debug
% s1_range = [1     :36000];
% s2_range = [36001 :48000];
% s3_range = [48001 :93600];
% s4_range = [93601 :129600];

%% filtering code

    highpass = 0.01; % 30 sec highpass filter
    Ca_fs    = 20; % 20hz
    
    ca_temp_data = ca_raw_data   ;      %　vlockedデータの代入 ca 縦：時間　ｘ　横：細胞
   
        [b,a] = butter(1, highpass /(Ca_fs/2),'high'); %
        %[b,a] = cheby2(1, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
        dataOut = filtfilt(b,a,ca_temp_data); %ゼロ位相filtfiltで ハイパス　フィルターを通す
        
%% zscoreも試す
	dataIn1 = zscore(dataOut(s1_range,:)); 	%　縦にフィルターがかかる
    dataIn2 = zscore(dataOut(s2_range,:)); 	%　縦にフィルターがかかる
    dataIn3 = zscore(dataOut(s3_range,:)); 	%　縦にフィルターがかかる
    dataIn4 = zscore(dataOut(s4_range,:)); 	%　縦にフィルターがかかる
    
        dataIn = [dataIn1; dataIn2; dataIn3; dataIn4];
        
%%
            ca_filt_data = (dataIn);

                negative = find(ca_filt_data<1); % 3SD以下をカット
    ca_filt_data(negative) = zeros(size(negative));
   
%% Save

          save (save_mat_name, 'ca_filt_data' )
              display('Finish saving');
    
end