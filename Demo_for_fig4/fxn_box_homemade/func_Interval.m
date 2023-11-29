function [peak_width] = func_Interval(ca_filt_temp,interval_range, thr_val) % ca_filt, interval_range = [1:1000], thr_val = cutoff criteria
%% comparison Interval        
%  ToIntervals(x,input) x=time stamp, input=logical 
 
% interval_range = [1:1000] ;  % input frame number
% thr_val = 3;

        A1 = ca_filt_temp(interval_range,:) > thr_val;
        time_stamp = [1:size(A1,1)]/20; % add time
    

        A1_table = [];
    
    for i_A1  = 1:size(A1,2)
        res_temp     = ToIntervals(time_stamp, A1(:,i_A1));
        A1_table = cat(1,A1_table,res_temp);
    end
    
    A1_table(:,3) = A1_table(:,2)-A1_table(:,1);
    zero_position = find(A1_table(:,3)==0);
    A1_table(zero_position,3) = 1/20; % 50ms‚Ì•‚Æ‚µ‚Ä0‚ğ’u‚«Š·‚¦
    peak_width = A1_table(:,3);    
end 
    
