%% MPPCA-ICA calculation code
clc; clear; close all; disp('Previous data is cleared'); tic;
data_cell = {}; shift=1; % Dont modify this line
addpath("fxn_box_homemade\"); 
load('mouse_ca_data.mat') % Ca data
%% load ca data 
ca_raw_data      = resv333py38thr1 ; % input time x neuron matrix

% Manuscript uses "bin_frame_num = 4", but it takas time hrs-half day. 
% To check code-behavior,  "bin_frame_num = 20" is recommended.

% bin_frame_num    = 4        ; % 20   -> 1s, 1s binning, For PCA-ICA, for Manuscript.
bin_frame_num    = 20        ; % 20   -> 1s, 1s binning, For PCA-ICA, I reccomend to not change here, first.
sample_fps       = 20        ; % RGECO: 10FPS, G-CaMP:20hz. Dependent on your data fps.
%%
data_cell{1+shift,1} = ('AB');   data_cell{1+shift,2} = [6498	:7699] ;   % REM_range_for_stat = [72 104; 379 604];
data_cell{2+shift,1} = ('BC');    data_cell{2+shift,2} = [7700    :8900] ;  
data_cell{3+shift,1} = ('CD');   data_cell{3+shift,2} = [8901   :10102] ; 
data_cell{4+shift,1} = ('DE');   data_cell{4+shift,2} = [10103	:11304] ; % REM_range_for_stat = [757 960];
data_cell{5+shift,1} = ('Rand');    data_cell{5+shift,2} = [11305   :16113] ; 
data_cell{6+shift,1} = ('T1');    data_cell{6+shift,2} = [16114   :17315] ; 
data_cell{7+shift,1}= ('awake');  data_cell{7+shift,2}= [17316	:18517] ; % [1000 1500]
data_cell{8+shift,1}= ('N2'); data_cell{8+shift,2}= [18518	:19718] ;
data_cell{9+shift,1}= ('R2'); data_cell{9+shift,2}= [19719	:20919] ;
data_cell{10+shift,1}= ('T2'); data_cell{10+shift,2}= [20920	:22120] ;
data_cell{11+shift,1}= ('T3'); data_cell{11+shift,2}= [22121	:23324] ;

%% data processing
% zscore: enable, cut-off: disable (Recommend zscore_mode:1, cut_off_mode:0)
filter_mode  = 0; % 0: disable, 1: enable 0.01hz-highpass butter worth filter.
zscore_mode  = 1; % 0: disable, 1: enable
cut_off_mode = 1; % 0: disable, 1: enable (If zscored, cut_off can be active.)
mu_mode      = 0; % 0: disable, 1: enable

[result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process_v3(data_cell, ... 
    ca_raw_data, bin_frame_num, sample_fps, filter_mode, zscore_mode, cut_off_mode, mu_mode);
%% # Assign reference and target section. Select either this section or below, (Disable eiteher here or below)
% # reference vs. all sesion analysis mode. I recomend this first to check data aspect.
reference_session_num = 10;
extra_session_num     = 10; % input 1st session num to exclude for backprojection. Latter session also will be removed. 
% target_session_num  = 1;

data_bin_reference = result_data_cell{reference_session_num+shift,7}; % input data for MPPCA-ICA
data_bin_target = result_ca_filt_data_ct; % backprojection to all data
% data_bin_target    = result_data_cell{target_session_num+shift,7}; % input data for MPPCA-ICA
%% Parameters
prms_MPPCA.prms_ca_fps                 = sample_fps;       % Confirm upper section you input. RGECO: 10FPS, G-CaMP:20hz. Dependent on your data fps.
prms_MPPCA.bin_frame_num               = bin_frame_num;    % Confirm upper section you input. 
prms_MPPCA.prms_forced_IC              = 0;        % def=0;    0:MP-PCA, 1-x:forced_IC num input;  
prms_MPPCA.prms_ICA_mode               = 3;        % def=3;    1:fastica, 2:fastICA, 3:Reconstruction ICA (RICA), 2 or 3 is better.
prms_MPPCA.prms_SD_thr                 = 2.5;      % def=2.5   (SD) ensemble weight threshold
prms_MPPCA.prms_reactivation_SD_thr    = 2;        % def=2;    (SD) filled coloring threshold for reactivation strength 
prms_MPPCA.prms_target_region1         = [501 : 1000 ] ;   % yellow, input second scale. like as [a:b]. disable: [0:0] 
prms_MPPCA.prms_target_region2         = [1501: 2000 ] ;   % green, input second scale. like as [c:d]. disable: [0:0]
prms_MPPCA.prms_ticklabel_mode         = 0 ;       % def=0;    0:off mode. 1:session name on mode. when you select whole data as target reference

% ### single mode for showing all figures ### Calculate MPPCA-ICA with backprojection onto all data
% [result_MPPCA] = fxn_MPPCA_ICAv08(data_bin_reference, data_bin_target, result_data_cell, prms_MPPCA); % calculation

%% (Double mode. Input two session values) Assign sessions

% input reference sessions
ref_session_num1 = 1; % select section
ref_session_num2 = 2; % select section
ref_session_num3 = 3; % select section
ref_session_num4 = 4; % select section

% input target session, which does not overlap with reference sessions
prms_multi.target_session_num  = 8 ; % Caution: Don't select same sassion to reference session, logically failed.
%% (Multi mode. Just run this section. Don't change) MP-PCA calculation
% 1st cal
data_bin_reference = result_data_cell{ref_session_num1+shift,7}; % input data for MPPCA-ICA
data_bin_target    = result_ca_filt_data_ct; % Don't change. for all-session backprojection.
[result_MPPCA1] = fxn_MPPCA_ICAv10(data_bin_reference, data_bin_target, result_data_cell, prms_MPPCA); % calculation
close all;% calculation

% 2nd cal
data_bin_reference = result_data_cell{ref_session_num2+shift,7}; % input data for MPPCA-ICA
data_bin_target    = result_ca_filt_data_ct; % Don't change. for all-session backprojection. 
[result_MPPCA2] = fxn_MPPCA_ICAv10(data_bin_reference, data_bin_target, result_data_cell, prms_MPPCA); % calculation
close all;% calculation

% 3rd cal
data_bin_reference = result_data_cell{ref_session_num3+shift,7}; % input data for MPPCA-ICA
data_bin_target    = result_ca_filt_data_ct; % Don't change. for all-session backprojection. 
[result_MPPCA3] = fxn_MPPCA_ICAv10(data_bin_reference, data_bin_target, result_data_cell, prms_MPPCA); % calculation
close all;% calculation

% 4th cal
data_bin_reference = result_data_cell{ref_session_num4+shift,7}; % input data for MPPCA-ICA
data_bin_target    = result_ca_filt_data_ct; % Don't change. for all-session backprojection. 
[result_MPPCA4] = fxn_MPPCA_ICAv10(data_bin_reference, data_bin_target, result_data_cell, prms_MPPCA); % calculation
close all;% calculation

% compose data
result_MPPCA_multi.result_MPPCA1 = result_MPPCA1;
result_MPPCA_multi.result_MPPCA2 = result_MPPCA2;
result_MPPCA_multi.result_MPPCA3 = result_MPPCA3;
result_MPPCA_multi.result_MPPCA4 = result_MPPCA4;
%% (Multi mode) Parameter input section
% Don't change
prms_multi.binarize_thr          = 2.5 ; % z-score, Calcium trace binalizing threshold.
prms_multi.shuffle_mode          = 2 ; % 1:randperm, 2:circshift
prms_multi.thr_percentile        = 1 ;  % def: 1%, top-bottom %-thresholding
prms_multi.iteration             = 1000;  % def: 1000, However, it take time 2-3 hrs.
prms_multi.ensemble_overlap_mode = 0;  % def:0, 1:allow to cal overlapping cells, but shuffling data can not function as control.
prms_multi.ensemble_set_A  = numel(result_MPPCA1.neuron_sig_IDs(:));
prms_multi.ensemble_set_B  = numel(result_MPPCA2.neuron_sig_IDs(:));
prms_multi.ensemble_set_C  = numel(result_MPPCA3.neuron_sig_IDs(:));
prms_multi.ensemble_set_D  = numel(result_MPPCA4.neuron_sig_IDs(:));
% prms_multi.target_session_num = is set upper section.
%% (Multi mode) Parameter input section
% ### input saved filename. To protect calculation resutls.

filename_for_save = (['Kareem_Results_of_connect_in_ref_1234_target5_cal','.mat']); 
%% (Multi mode) Start calculation. Need long time.
[result_final_cell, result_multi] = fxn_MPPCA_cal_connectivity_multi(result_data_cell, result_MPPCA_multi, prms_multi);
%%