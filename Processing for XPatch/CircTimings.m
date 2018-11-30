%%%
% File name: CircTimings.m
% Author: Calvin Kuo
% Date: 10/29/2018
%
% Pulls sensor sync timings. Designed to work with new firmware

function [offset, t, t1_begin, t2_begin, sensor1_interp, sensor2_interp] = CircTimings( sensor1_data, sensor2_data, startInd, startPage, nPages )
    
    sensor1_signal = [];
    sensor2_signal = [];
    
    curr_page = 0;
    currInd = startInd;
    currPage = startPage;
    totPage = 0;
    
    % Get Start Timing
    if currPage == 1
        t1_begin = sensor1_data(startInd).Info.t1 / 1000;
        t2_begin = sensor2_data(startInd).Info.t1 / 1000;
        currPage = 2;
        currSample = 22;
    elseif currPage == 2
        t1_begin = sensor1_data(startInd).Info.t2 / 1000;
        t2_begin = sensor2_data(startInd).Info.t2 / 1000;
        currPage = 3;
        currSample = 43;
    elseif currPage == 3
        t1_begin = sensor1_data(startInd).Info.t3 / 1000;
        t2_begin = sensor2_data(startInd).Info.t3 / 1000;
        currPage = 4;
        currSample = 64;
    else
        t1_begin = double(sensor1_data(startInd).Info.t4) / 1000;
        t2_begin = double(sensor2_data(startInd).Info.t4) / 1000;
        currPage = 1;
        currInd = currInd + 1;
        currSample = 1;
    end
    
    % Stitch signal together
    while( totPage < nPages )
        sensor1_signal = [sensor1_signal; ...
            sensor1_data(currInd).lin_acc(currSample,1) * 2 / 0.098 * 0.07 * pi / 180; ...
            sensor1_data(currInd).lin_acc(currSample,2) * 2 / 0.098 * 0.07 * pi / 180; ...
            sensor1_data(currInd).lin_acc(currSample,3) * 2 / 0.098 * 0.07 * pi / 180; ...
            sensor1_data(currInd).ang_vel(currSample,1); ...
            sensor1_data(currInd).ang_vel(currSample,2); ...
            sensor1_data(currInd).ang_vel(currSample,3)];
        
        sensor2_signal = [sensor2_signal; ...
            sensor2_data(currInd).lin_acc(currSample,1) * 2 / 0.098 * 0.07 * pi / 180; ...
            sensor2_data(currInd).lin_acc(currSample,2) * 2 / 0.098 * 0.07 * pi / 180; ...
            sensor2_data(currInd).lin_acc(currSample,3) * 2 / 0.098 * 0.07 * pi / 180; ...
            sensor2_data(currInd).ang_vel(currSample,1); ...
            sensor2_data(currInd).ang_vel(currSample,2); ...
            sensor2_data(currInd).ang_vel(currSample,3)];
        
        currSample = currSample + 1;
        if ( currSample == 22 )
            totPage = totPage + 1;
            currPage = currPage + 1;
        elseif ( currSample == 43 )
            totPage = totPage + 1;
            currPage = currPage + 1;
        elseif ( currSample == 64 )
            totPage = totPage + 1;
            currPage = currPage + 1;
        elseif ( currSample == 85 )
            totPage = totPage + 1;
            currPage = 1;
            currInd = currInd + 1;
            currSample = 1;
        end
    end
    
    % Get end timing
    if currPage == 1
        t1_end = double( sensor1_data(currInd-1).Info.t4 ) / 1000;
        t2_end = double( sensor2_data(currInd-1).Info.t4 ) / 1000;
    elseif currPage == 2
        t1_end = sensor1_data(currInd).Info.t1 / 1000;
        t2_end = sensor2_data(currInd).Info.t1 / 1000;
    elseif currPage == 3
        t1_end = sensor1_data(currInd).Info.t2 / 1000;
        t2_end = sensor2_data(currInd).Info.t2 / 1000;
    else
        t1_end = sensor1_data(currInd).Info.t3 / 1000;
        t2_end = sensor2_data(currInd).Info.t3 / 1000;
    end
    
    dt1 = ( t1_end - t1_begin ) / length( sensor1_signal );
    dt2 = ( t2_end - t2_begin ) / length( sensor2_signal );
    Fs1 = 1 / dt1;
    Fs2 = 1 / dt2;
    
    dt = 1/800;
    Fs = 800;
    
    t1 = ( 0:(length(sensor1_signal)-1) ) * dt1;
    t2 = ( 0:(length(sensor2_signal)-1) ) * dt2;
    
    t = 0:dt:min( [t1(end), t2(end)] );
    sensor1_interp = interp1( t1, sensor1_signal, t );
    sensor2_interp = interp1( t2, sensor2_signal, t );
    [acor, lag] = xcorr( sensor1_interp, sensor2_interp );
    [m,i] = max( acor );
    nOff = lag(i);
    if nOff < 0
        offset = t(-nOff) - t(1);
    elseif nOff > 0
        offset = t(1) - t(nOff);
    else
        offset = 0;
    end
end

function offset = CircTimingsObsolete( sensor1_data, sensor2_data, crange )

    t1_full = sensor1_data.t;
    t1init = find( t1_full >= crange(1), 1, 'first' );
    t1fin = find( t1_full >= crange(2), 1, 'first' );
    t1_cut = t1_full(t1init:t1fin);
    
    t2_full = sensor2_data.t;
    
    angvel1 = sqrt( sum( sensor1_data.ang_vel(t1init:t1fin,:) .* sensor1_data.ang_vel(t1init:t1fin,:), 2 ) );
    angvel2_raw = sqrt( sum( sensor2_data.ang_vel .* sensor2_data.ang_vel, 2 ) );
    angvel2 = interp1( t2_full, angvel2_raw, t1_cut );
    
    [acor, lag] = xcorr( angvel1, angvel2 );
    [m,i] = max( acor );
    nOff = lag(i);
    if nOff < 1
        offset = t1_cut(-nOff) - t1_cut(1);
    else
        offset = t1_cut(1) - t1_cut(nOff);
    end
end