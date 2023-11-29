%% MPPCA-ICA calculation code
clc; clear; close all; disp('Previous data is cleared'); tic;
data_cell = {}; global shift; shift=1; % Dont modify this line
addpath("fxn_box_homemade\"); 
load('TI-4') % Ca data
%% ### load ca data ###
ca_raw_data      = resv333py38thr1 ; % input time x neuron matrix
%% Input bahavioral session and extra session for reference.
data_cell{1+shift,1} = ('AB');   data_cell{1+shift,2} = [6498	:7699] ;   % REM_range_for_stat = [72 104; 379 604];
data_cell{2+shift,1} = ('BC');    data_cell{2+shift,2} = [7700    :8900] ;  
data_cell{3+shift,1} = ('CD');   data_cell{3+shift,2} = [8901   :10102] ; 
data_cell{4+shift,1} = ('DE');   data_cell{4+shift,2} = [10103	:11304] ; % REM_range_for_stat = [757 960];
data_cell{5+shift,1} = ('Rand');    data_cell{5+shift,2} = [11305   :16113] ; 
data_cell{6+shift,1} = ('T1');    data_cell{6+shift,2} = [16114   :17315] ; 
data_cell{7+shift,1}= ('awake');  data_cell{7+shift,2}= [17316	:18517] ; % [1000 1500]
data_cell{8+shift,1}= ('N2'); data_cell{8+shift,2}= [18518	:19718] ;
data_cell{9+shift,1}= ('R2'); data_cell{9+shift,2}= [19719	:20919] ;
data_cell{10+shift,1}= ('T2'); data_cell{8+shift,2}= [20920	:22120] ;
data_cell{11+shift,1}= ('T3'); data_cell{9+shift,2}= [22121	:23324] ;
%% data processing
% zscore: enable, cut-off: disable (Recommend zscore_mode:1, cut_off_mode:0)
prms_MPPCA.bin_frame_num  = 10  ; % 20   -> 1s, 1s binning, For PCA-ICA, I reccomend to not change here, first.
prms_MPPCA.sample_fps     = 20 ; % RGECO: 10FPS, G-CaMP:20hz. Dependent on your data fps.
prms_MPPCA.filter_mode    = 1  ; % 0: disable, 1: enable 0.01hz-highpass butter worth filter.
prms_MPPCA.zscore_mode    = 1  ; % 0: disable, 1: enable
prms_MPPCA.cut_off_mode   = 1  ; % 0: disable, 1: enable (If zscored, cut_off can be active.)
prms_MPPCA.mu_mode        = 0  ; % 0: disable, 1: enable

[result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process_v4(data_cell, ca_raw_data, prms_MPPCA);
%% Parameters
prms_MPPCA.prms_forced_IC              = 0;        % def=0;    0:MP-PCA, 1-x:forced_IC num input;  
prms_MPPCA.prms_ICA_mode               = 3;        % def=3;    1:fastica, 2:fastICA, 3:Reconstruction ICA (RICA), 2 or 3 is better.
prms_MPPCA.prms_SD_thr                 = 2.5;      % def=2.5   (SD) ensemble weight threshold
prms_MPPCA.prms_reactivation_SD_thr    = 2;        % def=2;    (SD) filled coloring threshold for reactivation strength 
prms_MPPCA.prms_target_region1         = [501 : 1000 ] ;   % new ver OFF, yellow, input second scale. like as [a:b]. disable: [0:0] 
prms_MPPCA.prms_target_region2         = [1501: 2000 ] ;   % new ver OFF, green, input second scale. like as [c:d]. disable: [0:0]
prms_MPPCA.prms_ticklabel_mode         = 1 ;       % def=0;    0:off mode. 1:session name on mode. when you select whole data as target reference
% Session info
prms_MPPCA.prms_reference_session_num  = 10 ; %
prms_MPPCA.prms_extra_session_mode     = 1 ; % 0:backpreoject to all session, 1:enable backprojection prior to extra session.
prms_MPPCA.prms_extra_session_num      = 10 ; % input 1st session num to exclude for backprojection. Latter session also will be removed. 
prms_MPPCA.prms_figure_all_mode        = 1 ; % 0:disable figure all, 1:enable figure all.

% ### single mode for showing all figures ### Calculate MPPCA-ICA with backprojection onto all data
[result_MPPCA] = fxn_MPPCA_ICAv11(result_ca_filt_data_ct, result_data_cell, prms_MPPCA); % calculation

%% Selective analysis and stat for backprojected reactivatin strength
data_ = result_MPPCA.r_strength_targetz; % for debug

% Input values with binned frame range for your interested behavioral session. 
select_cell{1+shift,1}= ('Selected-1'); select_cell{1+shift,2} = [6498	:7699] ; % REM_range_for_stat = [72 104; 379 604];
select_cell{2+shift,1}= ('Selected-2'); select_cell{2+shift,2} = [7700    :8900] ;  
select_cell{3+shift,1}= ('Selected-3'); select_cell{3+shift,2} = [8901   :10102] ; 
select_cell{4+shift,1}= ('Selected-4'); select_cell{4+shift,2} = [10103	:11304] ; % REM_range_for_stat = [757 960];
select_cell{5+shift,1}= ('Selected-5'); select_cell{5+shift,2} = [20920	:22120] ; 

% Input parameters
prms_MPPCA.prms_thr_freq_search = 2; % def=2 SD, SD threshold for frequency cal. McHugh uses Raw RS 5.
prms_MPPCA.prms_stat_ref        = 2; % 
prms_MPPCA.prms_stat_tar        = 3; % 
prms_MPPCA.prms_thr_stat_fold   = 1; % 
%% Wilcoxon ranksum stat for comparing patterns among original sessions
%  This code calculates by excluding extra session to avoid the change of entire back projected reactivation strength  
[result_MPPCA_original_stat] = fxn_MPPCA_selective_stat_original_frame(result_MPPCA, prms_MPPCA);

%% Wilcoxon ranksum stat for comparing patterns among selected and binned sessions 
%  This code calculates by excluding extra session to avoid the change of entire back projected reactivation strength  
[result_MPPCA_select_stat] = fxn_MPPCA_selective_stat_binned_frame(select_cell, result_MPPCA, prms_MPPCA);

%%