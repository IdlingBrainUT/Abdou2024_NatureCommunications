%% MPPCA-ICA calculation code
clc; clear; close all; disp('Previous data is cleared'); tic;
data_cell = {}; shift=1; % Dont modify this line
addpath("fxn_box_homemade\"); 
load('mouse_ca_data.mat') % Ca data

%% load ca data 
ca_raw_data      = resv333py38thr1; % input time x neuron matrix
bin_frame_num    = 4;   % 20   -> 1s, 1s binning, For PCA-ICA, I reccomend to not change here, first.

%% %%  mouse sod 5 session frame information

% input frame info
data_cell{1+shift,1} = ('Habit');   data_cell{1+shift,2} = [1	:4093] ;   % REM_range_for_stat = [72 104; 379 604];
data_cell{2+shift,1} = ('N1');    data_cell{2+shift,2} = [4094    :5295] ;  
data_cell{3+shift,1} = ('R1');   data_cell{3+shift,2} = [5296   :6497] ; 
data_cell{4+shift,1} = ('Training');   data_cell{4+shift,2} = [6498	:11304] ; % REM_range_for_stat = [757 960];
data_cell{5+shift,1} = ('Rand');    data_cell{5+shift,2} = [11305   :16113] ; 
data_cell{6+shift,1} = ('T1');    data_cell{6+shift,2} = [16114   :17315] ; 
data_cell{7+shift,1}= ('awake');  data_cell{7+shift,2}= [17316	:18517] ; % [1000 1500]
data_cell{8+shift,1}= ('N2'); data_cell{8+shift,2}= [18518	:19718] ;
data_cell{9+shift,1}= ('R2'); data_cell{9+shift,2}= [19719	:20919] ;
data_cell{10+shift,1}= ('T2'); data_cell{10+shift,2}= [20920	:22120] ;
data_cell{11+shift,1}= ('T3'); data_cell{11+shift,2}= [22121	:23324] ;

%% data processing
% This code includes z-score processing, even you select either "yes" or "no" to filter fluctuation.

% select "n".
[result_ca_filt_data_ct, result_data_cell] = fxn_MPPCA_ca_data_process(data_cell, ca_raw_data, bin_frame_num);

%% ### Following sessions are code for calculaitng the connectivity among different session patterns.
% ######  (Double mode)  ######
% Backprojection (Target session) should be assigned as all session
% To compare connectivity among patterns detected from differnet sessions, Time length of reactivation matrix must be identical. 
%% Parameters
prms_forced_IC              = 0;        % def=0;    0:MP-PCA, 1-x:forced_IC num input;  
prms_ICA_mode               = 3;        % def=3;    1:fastica, 2:fastICA, 3:Reconstruction ICA (RICA), 2 or 3 is better.
prms_SD_thr                 = 2;      % def=2.5   (SD) ensemble weight threshold
prms_reactivation_SD_thr    = 1.5;        % def=2;    (SD) filled coloring threshold for reactivation strength 
prms_target_region1         = [0:0] ;   % yellow, input second scale. like as [a:b]. disable: [0:0] 
prms_target_region2         = [0:0] ;   % green, input second scale. like as [c:d]. disable: [0:0]
prms_ticklabel_mode         = 1 ;       % def=0;    0:off mode. 1:session name on mode. when you select whole data as target reference
%% (Double mode. Input two session values) Assign sessions

reference_session_num1 = 6 ; % select section
reference_session_num2 = 8 ; % select section

prms.session_num_assess  = 9 ; % Caution: Don't select same sassion to reference session, logically failed.
%% (Double mode. Just run this section. Don't change) MP-PCA calculation
data_bin_z_reference = result_data_cell{reference_session_num1+shift,7}; % input data for MPPCA-ICA
data_bin_z_target    = result_ca_filt_data_ct; % Don't change. for all-session backprojection.

[result_MPPCA1] = fxn_MPPCA_ICAv05(data_bin_z_reference, data_bin_z_target, prms_forced_IC, prms_ICA_mode, prms_SD_thr,...
                                   prms_reactivation_SD_thr, prms_target_region1, prms_target_region2, bin_frame_num, ...
                                   result_data_cell, prms_ticklabel_mode); close all;% calculation

data_bin_z_reference = result_data_cell{reference_session_num2+shift,7}; % input data for MPPCA-ICA
data_bin_z_target    = result_ca_filt_data_ct; % Don't change. for all-session backprojection.

[result_MPPCA2] = fxn_MPPCA_ICAv05(data_bin_z_reference, data_bin_z_target, prms_forced_IC, prms_ICA_mode, prms_SD_thr,...
                                   prms_reactivation_SD_thr, prms_target_region1, prms_target_region2, bin_frame_num, ...
                                   result_data_cell, prms_ticklabel_mode); close all;% calculation%%                               

%% (Double mode) Assesing cell-cell pairwise interaction round-robin

% Don't change
prms.binarize_thr        = 3 ; % Calcium trace binalizing threshold.
prms.shuffle_mode        = 2 ;      % 1:randperm, 2:circshift
prms.ensemble_num1        = numel(result_MPPCA1.neuron_sig_IDs);
prms.ensemble_num2        = numel(result_MPPCA2.neuron_sig_IDs);
prms.thr_percentile      = 1; % def: 1%, top-bottom %-thresholding
%%
result_connectivity = fxn_MPPCA_ensemble_connectivity_double(result_data_cell, result_MPPCA1, result_MPPCA2, prms);
%% Figure PDF, matlab 2020a or later version required.
% Please wait for several minutes, because more than 10 by 10 figure takes time more than 5 min to draw figures even in high-spec PC.
extension = ('pdf'); % select: 'pdf' or 'jpg' for your normal figure, 'eps' for illustrator vector-format. 

fxn_MPPCA_figure_venn_histogram(result_connectivity, prms, extension); % figure Venn and Histogram
disp('Finish outputting figures. Now you can continue to use MATLAB') % Wait until this message appear.
%%
                               
                               