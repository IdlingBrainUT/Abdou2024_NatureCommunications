%% 
function [ca_digital] = fxn_mod_round_binning(ca_filt_data, bin_frame_num);
    mod_val = mod(size(ca_filt_data,1), bin_frame_num);
    
if mod_val == 0
    display('   mod_val = 0, normal binning done!')
    for i = 1:size(ca_filt_data,2)
    ca_mod(:,i)  = ca_filt_data(:,i);
    end
    [ca_mod_bin,  ca_mod_mean] = funcHF_temporal_binning(ca_mod, bin_frame_num);
else    
    display('      mod_val is not = 0, mod_rounded binning done! Last frame was rounded')
    for i = 1:size(ca_filt_data,2)
    ca_mod(:,i)  = [ca_filt_data(:,i); zeros(bin_frame_num-mod_val,1)];
    end 
    [ca_mod_bin,  ca_mod_mean] = funcHF_temporal_binning(ca_mod(1:end-bin_frame_num,:), bin_frame_num);
end
% [ca_mod_bin,  ca_mod_mean] = funcHF_temporal_binning(ca_mod, bin_frame_num);
    
    % for logical
%      ca_logical = ca_mod_bin > 0; 
%      ca_digital = double(ca_logical);  

% for num
ca_digital = ca_mod_bin; 

%      display('finish fxn_mod_binning')
end

%%