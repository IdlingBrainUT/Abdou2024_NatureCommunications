function ca_filt_data = funcHF_ca_filt_3z_f8(ca_raw_data, f1,f2,f3,f4,f5,f6,f7,f8);

%% filtering code

    highpass = 0.01; % 30 sec highpass filter
    Ca_fs    = 20; % 20hz
    
    ca_temp_data = ca_raw_data   ;      %　vlockedデータの代入 ca 縦：時間　ｘ　横：細胞
   
        [b,a] = butter(1, highpass /(Ca_fs/2),'high'); %
        %[b,a] = cheby2(1, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
        dataOut = filtfilt(b,a,ca_temp_data); %ゼロ位相filtfiltで ハイパス　フィルターを通す
        
%% zscoreも試す

% t.f1 = t.s1(1):t.s1(2);
% t.f2 = t.s2(1):t.s2(2);
% t.f3 = t.s3(1):t.s3(2);
% t.f4 = t.s4(1):t.s4_fr(2);
% t.f5 = t.s5(1):t.s5(2);
% t.f6 = t.s6(1):t.s6(2);
% t.f7 = t.s7(1):t.s7(2);
% t.f8 = t.s8(1):t.s8(2);

	dataIn1 = zscore(dataOut(f1,:)); 	    %　縦にフィルターがかかる
    dataIn2 = zscore(dataOut(f2,:)); 	%　縦にフィルターがかかる
    dataIn3 = zscore(dataOut(f3,:)); 	%　縦にフィルターがかかる
    dataIn4 = zscore(dataOut(f4,:)); 	%　縦にフィルターがかかる
    dataIn5 = zscore(dataOut(f5,:)); 	%　縦にフィルターがかかる
    dataIn6 = zscore(dataOut(f6,:)); 	%　縦にフィルターがかかる 
    dataIn7 = zscore(dataOut(f7,:)); 	%　縦にフィルターがかかる
    dataIn8 = zscore(dataOut(f8,:)); 	%　縦にフィルターがかかる 

        dataIn = [dataIn1; dataIn2; dataIn3; dataIn4; dataIn5; dataIn6; dataIn7; dataIn8];

           dataIn =  zscore(dataIn(:,:));
        
%%
            ca_filt_data = (dataIn);

                negative = find(ca_filt_data<2); % 2SD以下をカット
    ca_filt_data(negative) = zeros(size(negative));
    
%         % 3SD後に -3 でベースラインに戻すコード
%     ca_filt_data_zero = ca_filt_data-3;                
%     negative = find(ca_filt_data_zero<3);
%     ca_filt_data_zero(negative) = zeros(size(negative));
   
%%
    
end