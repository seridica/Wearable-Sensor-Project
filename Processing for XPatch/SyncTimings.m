%%%
% File name: SyncTimings.m
% Author: Calvin Kuo
% Date: 10/29/2018
%
% Pulls sensor sync timings

function timings = SyncTimings( sensor_data, num_syncs, start_ind, offset )
    this_ind = start_ind;
    
    timings = zeros( num_syncs, 1 );
    
    for i=1:num_syncs
        
        t1 = sensor_data(this_ind).Info.t1;
        t3 = sensor_data(this_ind).Info.t3;
        
        if (offset == 0)
            offset = offset + 1;
            timings(i) = t1;
        else
            offset = 0;
            this_ind = this_ind + 1;
            timings(i) = t3;
        end
    end
end