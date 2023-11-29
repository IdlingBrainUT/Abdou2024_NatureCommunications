function ca_filt_data = fxn_ca_HPfilt_15s(ca_raw_data,t);

%% filtering code

    highpass = 0.01; % 30 sec highpass filter
    Ca_fs    = 20; % 20hz
    
    ca_temp_data = ca_raw_data   ;      %　vlockedデータの代入 ca 縦：時間　ｘ　横：細胞
   
        [b,a] = butter(1, highpass /(Ca_fs/2),'high'); %
        %[b,a] = cheby2(1, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
%         dataOut = filtfilt(b,a,ca_temp_data); %ゼロ位相filtfiltで ハイパス　フィルターを通す
        dataOut = ca_temp_data;
%% zscoreも試す
	dataIn1 = zscore(dataOut(t.s1,:)); 	    %　縦にフィルターがかかる
    dataIn2 = zscore(dataOut(t.s2,:)); 	%　縦にフィルターがかかる
    dataIn3 = zscore(dataOut(t.s3,:)); 	%　縦にフィルターがかかる
    dataIn4 = zscore(dataOut(t.s4,:)); 	%　縦にフィルターがかかる
    dataIn5 = zscore(dataOut(t.s5,:)); 	%　縦にフィルターがかかる
    dataIn6 = zscore(dataOut(t.s6,:)); 	%　縦にフィルターがかかる
    dataIn7 = zscore(dataOut(t.s7,:)); 	%　縦にフィルターがかかる
    dataIn8 = zscore(dataOut(t.s8,:)); 	%　縦にフィルターがかかる
    dataIn9 = zscore(dataOut(t.s9,:)); 	%　縦にフィルターがかかる
    dataIn10 = zscore(dataOut(t.s10,:)); 	%　縦にフィルターがかかる
    dataIn11 = zscore(dataOut(t.s11,:)); 	%　縦にフィルターがかかる
    dataIn12 = zscore(dataOut(t.s12,:)); 	%　縦にフィルターがかかる
    dataIn13 = zscore(dataOut(t.s13,:)); 	%　縦にフィルターがかかる
    dataIn14 = zscore(dataOut(t.s14,:)); 	%　縦にフィルターがかかる
    dataIn15 = zscore(dataOut(t.s15,:)); 	%　縦にフィルターがかかる
 
    
        dataIn = [dataIn1; dataIn2; dataIn3; dataIn4; dataIn5; dataIn6; ...
                  dataIn7; dataIn8; dataIn9; dataIn10; dataIn11; dataIn12; dataIn13; dataIn14; dataIn15];

%            dataIn =  zscore(dataOut(:,:));
        
%%
            ca_filt_data = (dataIn);

                negative = find(ca_filt_data<3); % 2SD以下をカット
    ca_filt_data(negative) = zeros(size(negative));
    
%         % 3SD後に -3 でベースラインに戻すコード
%     ca_filt_data_zero = ca_filt_data-3;                
%     negative = find(ca_filt_data_zero<3);
%     ca_filt_data_zero(negative) = zeros(size(negative));
   
%%
    
end