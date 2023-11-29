function ca_filt_data = fxn_ca_HPfilt_15s(ca_raw_data,t);

%% filtering code

    highpass = 0.01; % 30 sec highpass filter
    Ca_fs    = 20; % 20hz
    
    ca_temp_data = ca_raw_data   ;      %�@vlocked�f�[�^�̑�� ca �c�F���ԁ@���@���F�זE
   
        [b,a] = butter(1, highpass /(Ca_fs/2),'high'); %
        %[b,a] = cheby2(1, 10, highpass /(Ca_fs/2), 'high');% highpass =0.01 good  
%         dataOut = filtfilt(b,a,ca_temp_data); %�[���ʑ�filtfilt�� �n�C�p�X�@�t�B���^�[��ʂ�
        dataOut = ca_temp_data;
%% zscore������
	dataIn1 = zscore(dataOut(t.s1,:)); 	    %�@�c�Ƀt�B���^�[��������
    dataIn2 = zscore(dataOut(t.s2,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn3 = zscore(dataOut(t.s3,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn4 = zscore(dataOut(t.s4,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn5 = zscore(dataOut(t.s5,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn6 = zscore(dataOut(t.s6,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn7 = zscore(dataOut(t.s7,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn8 = zscore(dataOut(t.s8,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn9 = zscore(dataOut(t.s9,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn10 = zscore(dataOut(t.s10,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn11 = zscore(dataOut(t.s11,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn12 = zscore(dataOut(t.s12,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn13 = zscore(dataOut(t.s13,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn14 = zscore(dataOut(t.s14,:)); 	%�@�c�Ƀt�B���^�[��������
    dataIn15 = zscore(dataOut(t.s15,:)); 	%�@�c�Ƀt�B���^�[��������
 
    
        dataIn = [dataIn1; dataIn2; dataIn3; dataIn4; dataIn5; dataIn6; ...
                  dataIn7; dataIn8; dataIn9; dataIn10; dataIn11; dataIn12; dataIn13; dataIn14; dataIn15];

%            dataIn =  zscore(dataOut(:,:));
        
%%
            ca_filt_data = (dataIn);

                negative = find(ca_filt_data<3); % 2SD�ȉ����J�b�g
    ca_filt_data(negative) = zeros(size(negative));
    
%         % 3SD��� -3 �Ńx�[�X���C���ɖ߂��R�[�h
%     ca_filt_data_zero = ca_filt_data-3;                
%     negative = find(ca_filt_data_zero<3);
%     ca_filt_data_zero(negative) = zeros(size(negative));
   
%%
    
end