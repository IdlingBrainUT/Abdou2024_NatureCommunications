function [filt_raw_lfp] = fxn_bandstop_filt_lfp(raw_lfp_data, fs, fpass)
%% for debug,
% y = bandstop(x,fpass,fs); % x: signal, fpass = [hz hz], fs
% raw_lfp_data = raw_LFP1;
% fpass = [40 80]; % 60 +- 15 hz
% fs = fs_LFP;
%% filtering code
filt_raw_lfp = bandstop(raw_lfp_data,fpass,fs,'Steepness',0.95);
%% figure
% figure;
% h1 = subplot(311); plot(raw_lfp_data,'b'); 
% h2 = subplot(312); plot(filt_raw_lfp,'r'); 
% h3 = subplot(313); plot(raw_lfp_data,'b'); hold on;
%               plot(filt_raw_lfp,'r'); 
%               
% linkaxes([h1,h2,h3],'x' )
%%    
end