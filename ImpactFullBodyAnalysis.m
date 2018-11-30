%%%
% File name: ImpactFullBodyAnalysis.m
% Author: Calvin Kuo
% Date: 10/29/2018
%
% Code for performing a more full body analysis on impacts.

function [tplt, impact_linacc, impact_angvel, impact_angacc] = ImpactFullBodyAnalysis( data_inputs, rotMatrices, timing_offsets, t_start, t_end, linacc_thresh, wrange, imppre, imppost )

    % Close plots
    close all;
    
    nInputs = length( data_inputs );
    interpFs = 400;
    interpdt = 1/interpFs;
    
    %% Unified time vector
    tInit = -inf;
    tFin = inf;
    for i=1:length( data_inputs )
        curr_data = data_inputs(i);
        if ( curr_data.t(1) > tInit )
            tInit = curr_data.t(1);
        end
        if ( curr_data.t(end) < tFin )
            tFin = curr_data.t(end);
        end
    end
    
    %% Find Impacts
    % Use the first sensor data as the timing reference
    input_data = data_inputs(1);
    t_full = input_data.t;
    linacc_fullraw = input_data.lin_acc;
    
    gyroDelay = 0.0;
    t_base = tInit:interpdt:(tFin-gyroDelay); %t_full(1):interpdt:t_full(end);
    linacc_interp = interp1( t_full, linacc_fullraw, t_base, 'spline', 'extrap' );

    % Filter
    [b,a] = butter( 2, 190 / (interpFs/2) );
    linacc_interpfilt = filtfilt( b, a, linacc_interp );
    
    % Find event indices
    event_inds = FindPeaks( t_base, linacc_interpfilt, t_start, t_end, linacc_thresh, wrange );
    unique_inds = unique( event_inds );
    
    %% Cycle through data inputs and plot all on single plots
    nImpacts = length( unique_inds );
    
    % Limit to 50 for now
    if nImpacts > 5
        nImpacts = 5;
    end
    
    % Cycle data
    tplt = ( -imppre:imppost ) * interpdt;
    
    for i=1:nInputs
        
        % Pull full data and process
        input_data = data_inputs(i);
        t_full = input_data.t - timing_offsets(i);
        linacc_fullraw = input_data.lin_acc;
        angvel_fullraw = input_data.ang_vel;
        
        sensRot = rotMatrices(:,:,i);
        linacc_fullrot = zeros( size( linacc_fullraw ) );
        angvel_fullrot = zeros( size( angvel_fullraw ) );
        for j=1:length( t_full )
            linacc_fullrot(j,:) = ( sensRot * [-1 0 0; 0 -1 0; 0 0 1] * linacc_fullraw(j,:)' )';  % lin accel sensor frame to gyro sensor frame
            angvel_fullrot(j,:) = ( sensRot * angvel_fullraw(j,:)' )';
        end
        
        % Timing based on t
        dt = ( t_full(end) - t_full(1) ) / length( t_full );
        i
        Fs = 1/dt
        
        linacc_interp = interp1( t_full, linacc_fullrot, t_base, 'spline', 'extrap' );
        angvel_interp = interp1( t_full-gyroDelay, angvel_fullrot, t_base, 'spline', 'extrap' );
        
        % Filter
        [b,a] = butter( 2, 190 / (interpFs/2) );
        angvel_interpfilt = filtfilt( b, a, angvel_interp );
        linacc_interpfilt = filtfilt( b, a, linacc_interp );
        angacc_interpfilt = dt_order_4_or_five_point_stencil( t_base, angvel_interpfilt' )';
        
        % Pull impacts
        figure(100); % Temporary figure for plotting
        linacc_impacts = HelpPlotEventTraces( interpdt, linacc_interpfilt, unique_inds, imppre, imppost );
        angvel_impacts = HelpPlotEventTraces( interpdt, angvel_interpfilt, unique_inds, imppre, imppost );
        angacc_impacts = HelpPlotEventTraces( interpdt, angacc_interpfilt, unique_inds, imppre, imppost );
        
        impact_linacc(i) = linacc_impacts;
        impact_angvel(i) = angvel_impacts;
        impact_angacc(i) = angacc_impacts;
        
        clf;
        for j=1:nImpacts
            figure(j);
            subplot(2,2,1); hold on;
            plot( tplt, linacc_impacts.signal_x_impacts(:,j) );
            subplot(2,2,2); hold on;
            plot( tplt, linacc_impacts.signal_y_impacts(:,j) );
            subplot(2,2,3); hold on;
            plot( tplt, linacc_impacts.signal_z_impacts(:,j) );
            subplot(2,2,4); hold on;
            plot( tplt, linacc_impacts.signal_mag_impacts(:,j) );
            
            figure(j+nImpacts+1);
            subplot(2,2,1); hold on;
            plot( tplt, angvel_impacts.signal_x_impacts(:,j) );
            subplot(2,2,2); hold on;
            plot( tplt, angvel_impacts.signal_y_impacts(:,j) );
            subplot(2,2,3); hold on;
            plot( tplt, angvel_impacts.signal_z_impacts(:,j) );
            subplot(2,2,4); hold on;
            plot( tplt, angvel_impacts.signal_mag_impacts(:,j) );
            
            figure(j+2*nImpacts+2);
            subplot(2,2,1); hold on;
            plot( tplt, angacc_impacts.signal_x_impacts(:,j) );
            subplot(2,2,2); hold on;
            plot( tplt, angacc_impacts.signal_y_impacts(:,j) );
            subplot(2,2,3); hold on;
            plot( tplt, angacc_impacts.signal_z_impacts(:,j) );
            subplot(2,2,4); hold on;
            plot( tplt, angacc_impacts.signal_mag_impacts(:,j) );
        end
        figure( nImpacts+1 );
        subplot(2,2,1); hold on;
        plot( tplt, mean( linacc_impacts.signal_x_impacts, 2 ) );
        subplot(2,2,2); hold on;
        plot( tplt, mean( linacc_impacts.signal_y_impacts, 2 ) );
        subplot(2,2,3); hold on;
        plot( tplt, mean( linacc_impacts.signal_z_impacts, 2 ) );
        subplot(2,2,4); hold on;
        plot( tplt, mean( linacc_impacts.signal_mag_impacts, 2 ) );
        
        figure( nImpacts*2+2 );
        subplot(2,2,1); hold on;
        plot( tplt, mean( angvel_impacts.signal_x_impacts, 2 ) );
        subplot(2,2,2); hold on;
        plot( tplt, mean( angvel_impacts.signal_y_impacts, 2 ) );
        subplot(2,2,3); hold on;
        plot( tplt, mean( angvel_impacts.signal_z_impacts, 2 ) );
        subplot(2,2,4); hold on;
        plot( tplt, mean( angvel_impacts.signal_mag_impacts, 2 ) );
        
        figure(3*nImpacts+3);
        subplot(2,2,1); hold on;
        plot( tplt, mean( angacc_impacts.signal_x_impacts, 2 ) );
        subplot(2,2,2); hold on;
        plot( tplt, mean( angacc_impacts.signal_y_impacts, 2 ) );
        subplot(2,2,3); hold on;
        plot( tplt, mean( angacc_impacts.signal_z_impacts, 2 ) );
        subplot(2,2,4); hold on;
        plot( tplt, mean( angacc_impacts.signal_mag_impacts, 2 ) );
    end
end