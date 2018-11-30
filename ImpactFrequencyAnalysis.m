%%%
% File name: ImpactFrequencyAnalysis.m
% Author: Calvin Kuo
% Date: 10/19/2018
%
% Code for doing frequency analysis on IMU data from various events of
% interest. Contains functions for a number of separate analyses.

function ImpactFrequencyAnalysis( input_data, t_start, t_end, linacc_thresh, wrange, imppre, imppost, lin_ind, ang_ind )

    % Close plots
    close all;

    %r = [0.05; -0.03; -0.08];
    r = [0; 0; 0;];
    
    % Full data pull
    t_full = input_data.t;
    linacc_fullraw = input_data.lin_acc;
    angvel_fullraw = input_data.ang_vel;

    % Timing based on t
    dt = ( t_full(end) - t_full(1) ) / length( t_full );
    Fs = 1/dt;
    
    % Filter
    [b,a] = butter( 2, 190 / (Fs/2) );
    angvel_fullfilt = filtfilt( b, a, angvel_fullraw );
    angacc_fullfilt = dt_order_4_or_five_point_stencil( t_full, angvel_fullfilt' )';
    linacc_fullraw = filtfilt( b, a, linacc_fullraw );
    linacc_fullfilt = zeros( size( linacc_fullraw ) );
    for j=1:length( t_full )
        linacc_fullfilt(j,:) = ( linacc_fullraw(j,:)'*9.81 + cross( angacc_fullfilt(j,:)', r ) + cross( angvel_fullfilt(j,:)', cross( angvel_fullfilt(j,:)', r ) ) )' / 9.81;
    end
    
    % Find event indices
    %event_inds = FindPeaks( t_full, linacc_fullfilt, t_start, t_end, linacc_thresh, wrange );
    if lin_ind ~= 0
        event_inds = FindPeaks( t_full, -linacc_fullfilt(:,lin_ind), t_start, t_end, linacc_thresh, wrange );
    else
        event_inds = FindPeaks( t_full, linacc_fullfilt, t_start, t_end, linacc_thresh, wrange );
    end
    unique_inds = unique( event_inds );
    
    % Plots of traces
    figure(1); clf; hold on;
    linacc_impacts = HelpPlotEventTraces(dt, linacc_fullfilt, unique_inds, imppre, imppost);
    
    figure(2); clf; hold on;
    angvel_impacts = HelpPlotEventTraces(dt, angvel_fullfilt, unique_inds, imppre, imppost);
    
    figure(3); clf; hold on;
    angacc_impacts = HelpPlotEventTraces(dt, angacc_fullfilt, unique_inds, imppre, imppost);
    
    % Wavelet transforms
    figure(4); clf; hold on;
    EventAllWavelets( Fs, linacc_impacts );
    
    figure(5); clf; hold on;
    EventAllWavelets( Fs, angvel_impacts );
    
    figure(6); clf; hold on;
    EventAllWavelets( Fs, angacc_impacts );
    
    % Attenuation
    freqs = [150, 100, 60, 30, 10. 5];
    figure(7); clf; hold on;
    AttenuationPlot( linacc_fullfilt(:,lin_ind), Fs, freqs, unique_inds, imppre, imppost );
    
    figure(8); clf; hold on;
    AttenuationPlot( angvel_fullfilt(:,ang_ind), Fs, freqs, unique_inds, imppre, imppost );
    
    figure(9); clf; hold on;
    AttenuationPlot( angacc_fullfilt(:,ang_ind), Fs, freqs, unique_inds, imppre, imppost );
    
    % Normal FFT
    figure(10); clf; hold on;
    HelperPlotFft( t_full, linacc_fullfilt(:,lin_ind), t_start, t_end, 150 );
    
    figure(11); clf; hold on;
    HelperPlotFft( t_full, angvel_fullfilt(:,ang_ind), t_start, t_end, 150 );
    
    figure(12); clf; hold on;
    HelperPlotFft( t_full, angacc_fullfilt(:,ang_ind), t_start, t_end, 150 );
    %HelpLaurens( dt, angvel_fullfilt, angacc_fullfilt, linacc_fullfilt, unique_inds, imppre, imppost );
    %HelpLaurensFull( t_full, angvel_fullfilt, angacc_fullfilt, linacc_fullfilt, t_start, t_end );
end

%% Helper function for wavelet transforms of the impacts
function EventAllWavelets( Fs, event_signals )
    subplot(2,2,4); hold on;
    peak_percent = max( mean( event_signals.signal_mag_impacts' ) );
    maxc = EventWavelets( Fs, event_signals.signal_mag_impacts, peak_percent );
    xlabel( 'Time (s)' )
    ylabel( 'Frequency (Hz)' )
    set(gca,'FontSize',20)
    
    subplot(2,2,1); hold on;
    EventWavelets( Fs, event_signals.signal_x_impacts, peak_percent, [0, maxc] );
    xlabel( 'Time (s)' )
    ylabel( 'Frequency (Hz)' )
    set(gca,'FontSize',20)
    
    subplot(2,2,2); hold on;
    EventWavelets( Fs, event_signals.signal_y_impacts, peak_percent, [0, maxc] );
    xlabel( 'Time (s)' )
    ylabel( 'Frequency (Hz)' )
    set(gca,'FontSize',20)
    
    subplot(2,2,3); hold on;
    EventWavelets( Fs, event_signals.signal_z_impacts, peak_percent, [0, maxc] );
    xlabel( 'Time (s)' )
    ylabel( 'Frequency (Hz)' )
    set(gca,'FontSize',20)
end

%% Helper function for Laurens Model run
function HelpLaurens( dt, angvel_full, angacc_full, linacc_full, unique_inds, imppre, imppost )
    [b,a] = butter(2, 10 / (1/dt/2) );
    angvel_filt = filtfilt(b,a,angvel_full);
    angacc_filt = filtfilt(b,a,angacc_full);
    linacc_filt = filtfilt(b,a,linacc_full);
    t = ( 0:(imppre + imppost) )' * dt;
    
    % Filtered
    angvel_x_impacts = PullImpacts( angvel_filt(:,1), unique_inds, imppre, imppost );
    angvel_y_impacts = PullImpacts( angvel_filt(:,2), unique_inds, imppre, imppost );
    angvel_z_impacts = PullImpacts( angvel_filt(:,3), unique_inds, imppre, imppost );
    angvel_filtered = [ mean( angvel_x_impacts' )', mean( angvel_y_impacts' )', mean( angvel_z_impacts' )' ];

    angacc_x_impacts = PullImpacts( angacc_filt(:,1), unique_inds, imppre, imppost );
    angacc_y_impacts = PullImpacts( angacc_filt(:,2), unique_inds, imppre, imppost );
    angacc_z_impacts = PullImpacts( angacc_filt(:,3), unique_inds, imppre, imppost );
    angacc_filtered = [ mean( angacc_x_impacts' )', mean( angacc_y_impacts' )', mean( angacc_z_impacts' )' ];
    %angacc_filtered = dt_order_4_or_five_point_stencil( t', angvel_filtered' )';
    
    linacc_x_impacts = PullImpacts( linacc_filt(:,1), unique_inds, imppre, imppost );
    linacc_y_impacts = PullImpacts( linacc_filt(:,2), unique_inds, imppre, imppost );
    linacc_z_impacts = PullImpacts( linacc_filt(:,3), unique_inds, imppre, imppost );
    linacc_filtered = [ mean( linacc_x_impacts' )', mean( linacc_y_impacts' )', mean( linacc_z_impacts' )' ];
    
    t = ( 0:(imppre + imppost) )' * dt;
    zz = zeros( length(t), 1 );
    
    % Angular velocity input
    angvel_x_impacts = PullImpacts( angvel_full(:,1), unique_inds, imppre, imppost );
    angvel_y_impacts = PullImpacts( angvel_full(:,2), unique_inds, imppre, imppost );
    angvel_z_impacts = PullImpacts( angvel_full(:,3), unique_inds, imppre, imppost );
    angvel_impacts = [ mean( angvel_x_impacts' )', mean( angvel_y_impacts' )', mean( angvel_z_impacts' )' ];
    %angvel_impacts = [zz, zz, zz];
    
    % Angular Acceleration input
    angacc_x_impacts = PullImpacts( angacc_full(:,1), unique_inds, imppre, imppost );
    angacc_y_impacts = PullImpacts( angacc_full(:,2), unique_inds, imppre, imppost );
    angacc_z_impacts = PullImpacts( angacc_full(:,3), unique_inds, imppre, imppost );
    angacc_impacts = [ mean( angacc_x_impacts' )', mean( angacc_y_impacts' )', mean( angacc_z_impacts' )' ];
    %angacc_impacts = dt_order_4_or_five_point_stencil( t', angvel_impacts' )';
    
    % Linear Acceleration input
    linacc_x_impacts = PullImpacts( linacc_full(:,1), unique_inds, imppre, imppost );
    linacc_y_impacts = PullImpacts( linacc_full(:,2), unique_inds, imppre, imppost );
    linacc_z_impacts = PullImpacts( linacc_full(:,3), unique_inds, imppre, imppost );
    linacc_impacts = [ mean( linacc_x_impacts' )', mean( linacc_y_impacts' )', mean( linacc_z_impacts' )' ];
    %linacc_impacts = [zz,zz,zz];
    
    [ve_full, vsf_full] = RunLaurens( t, angvel_impacts, angacc_impacts, linacc_impacts );
    [ve_filt, vsf_filt] = RunLaurens( t, angvel_filtered, angacc_filtered, linacc_filtered );
    subplot(3,3,1); hold on;
    plot( t, angvel_impacts(:,1), 'k', 'LineWidth', 2.0 );
    plot( ve_full(:,1), ve_full(:,2) );
    plot( ve_filt(:,1), ve_filt(:,2) );
    legend( 'Measured', 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,2); hold on;
    plot( t, angvel_impacts(:,1) - ve_full(:,2) );
    plot( t, angvel_impacts(:,1) - ve_filt(:,2) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,3); hold on;
    plot( t, vsf_full(:,2) );
    plot( t, vsf_filt(:,2) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,4); hold on;
    plot( t, angvel_impacts(:,2), 'k', 'LineWidth', 2.0 );
    plot( ve_full(:,1), ve_full(:,3) );
    plot( ve_filt(:,1), ve_filt(:,3) );
    legend( 'Measured', 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,5); hold on;
    plot( ve_full(:,1), angvel_impacts(:,2) - ve_full(:,3) );
    plot( ve_filt(:,1), angvel_impacts(:,2) - ve_filt(:,3) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,6); hold on;
    plot( t, vsf_full(:,3) );
    plot( t, vsf_filt(:,3) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,7); hold on;
    plot( t, angvel_impacts(:,3), 'k', 'LineWidth', 2.0 );
    plot( ve_full(:,1), ve_full(:,4) );
    plot( ve_filt(:,1), ve_filt(:,4) );
    legend( 'Measured', 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,8); hold on;
    plot( ve_full(:,1), angvel_impacts(:,3) - ve_full(:,4) );
    plot( ve_filt(:,1), angvel_impacts(:,3) - ve_filt(:,4) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,9); hold on;
    plot( t, vsf_full(:,4) );
    plot( t, vsf_filt(:,4) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
end

%% Helper function for Laurens Model Run FUll
function HelpLaurensFull( t, angvel_full, angacc_full, linacc_full, t_start, t_end )
    % Pull section of signal and time
    tinit = find( t > t_start, 1, 'first' );
    tfin = find( t > t_end, 1, 'first' );
    
    [b,a] = butter(2, 10 / 200);
    angvel_filt = filtfilt(b,a,angvel_full);
    angacc_filt = filtfilt(b,a,angacc_full);
    linacc_filt = filtfilt(b,a,linacc_full);
    
    % Make magnitude if not already magnitude
    angvel_ss = angvel_full( tinit:tfin, : );
    angacc_ss = angacc_full( tinit:tfin, : );
    linacc_ss = linacc_full( tinit:tfin, : );
    
    angvel_ssf = angvel_filt( tinit:tfin, : );
    angacc_ssf = angacc_filt( tinit:tfin, : );
    linacc_ssf = linacc_filt( tinit:tfin, : );
    
    [ve_full, vsf_full] = RunLaurens( t(tinit:tfin)', angvel_ss, angacc_ss, linacc_ss );
    [ve_filt, vsf_filt] = RunLaurens( t(tinit:tfin)', angvel_ssf, angacc_ssf, linacc_ssf );
    subplot(3,3,1); hold on;
    plot( t(tinit:tfin), angvel_ss(:,1), 'k', 'LineWidth', 2.0 );
    plot( ve_full(:,1), ve_full(:,2) );
    plot( ve_filt(:,1), ve_filt(:,2) );
    legend( 'Measured', 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,2); hold on;
    plot( t(tinit:tfin), angvel_ss(:,1) - ve_full(:,2) );
    plot( t(tinit:tfin), angvel_ss(:,1) - ve_filt(:,2) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,3); hold on;
    plot( t(tinit:tfin), vsf_full(:,2) );
    plot( t(tinit:tfin), vsf_filt(:,2) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,4); hold on;
    plot( t(tinit:tfin), angvel_ss(:,2), 'k', 'LineWidth', 2.0 );
    plot( ve_full(:,1), ve_full(:,3) );
    plot( ve_filt(:,1), ve_filt(:,3) );
    legend( 'Measured', 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,5); hold on;
    plot( ve_full(:,1), angvel_ss(:,2) - ve_full(:,3) );
    plot( ve_filt(:,1), angvel_ss(:,2) - ve_filt(:,3) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,6); hold on;
    plot( t(tinit:tfin), vsf_full(:,3) );
    plot( t(tinit:tfin), vsf_filt(:,3) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,7); hold on;
    plot( t(tinit:tfin), angvel_ss(:,3), 'k', 'LineWidth', 2.0 );
    plot( ve_full(:,1), ve_full(:,4) );
    plot( ve_filt(:,1), ve_filt(:,4) );
    legend( 'Measured', 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,8); hold on;
    plot( ve_full(:,1), angvel_ss(:,3) - ve_full(:,4) );
    plot( ve_filt(:,1), angvel_ss(:,3) - ve_filt(:,4) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    subplot(3,3,9); hold on;
    plot( t(tinit:tfin), vsf_full(:,4) );
    plot( t(tinit:tfin), vsf_filt(:,4) );
    legend( 'Unfiltered Input', '10Hz Filtered Input' )
    
    total_error_full_abs = [ sum( abs( angvel_ss(:,1) - ve_full(:,2) ) ), ...
                    sum( abs( angvel_ss(:,2) - ve_full(:,3) ) ), ...
                    sum( abs( angvel_ss(:,3) - ve_full(:,4) ) ) ] / length( t(tinit:tfin) )
                
    total_error_filt_abs = [ sum( abs( angvel_ss(:,1) - ve_filt(:,2) ) ), ...
                    sum( abs( angvel_ss(:,2) - ve_filt(:,3) ) ), ...
                    sum( abs( angvel_ss(:,3) - ve_filt(:,4) ) ) ] / length( t(tinit:tfin) )
                
    total_error_full = [ sum( angvel_ss(:,1) - ve_full(:,2) ), ...
                    sum( angvel_ss(:,2) - ve_full(:,3) ), ...
                    sum( angvel_ss(:,3) - ve_full(:,4) ) ] / length( t(tinit:tfin) )
                
    total_error_filt = [ sum( angvel_ss(:,1) - ve_filt(:,2) ), ...
                    sum( angvel_ss(:,2) - ve_filt(:,3) ), ...
                    sum( angvel_ss(:,3) - ve_filt(:,4) ) ] / length( t(tinit:tfin) )
                
    total_storage_full = [ sum( abs( vsf_full(:,2) ) ), ...
                           sum( abs( vsf_full(:,3) ) ), ...
                           sum( abs( vsf_full(:,4) ) ) ] / length( t(tinit:tfin) )
                       
    total_storage_filt = [ sum( abs( vsf_filt(:,2) ) ), ...
                           sum( abs( vsf_filt(:,3) ) ), ...
                           sum( abs( vsf_filt(:,4) ) ) ] / length( t(tinit:tfin) )
end