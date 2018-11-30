%%%
% File name: AltSyncTimings.m
% Author: Calvin Kuo
% Date: 10/29/2018
%
% Pulls sensor sync timings

function timings = AltSyncTimings( sensor_data, crange )

    t_full = sensor_data.t;
    tinit = find( t_full >= crange(1), 1, 'first' );
    tfin = find( t_full >= crange(2), 1, 'first' );
    t_cut = t_full(tinit:tfin);
    
    angvel = sqrt( sum( sensor_data.ang_vel(tinit:tfin,:) .* sensor_data.ang_vel(tinit:tfin,:), 2 ) );
    [pks, locs] = findpeaks(angvel,'MinPeakDistance',100, 'MinPeakProminence',4);
    timings = t_cut(locs);
end