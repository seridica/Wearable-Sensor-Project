%% Helper plot function for impact time series
function func_out = HelpPlotEventTraces( dt, signal_full, event_inds, imppre, imppost, fignum )
    
    % If the signal you pass in is a triaxial signal, then creates 4
    % subplots with each axis and magnitude. If you pass a single signal,
    % then just plots that signal.
    t = ( 0:(imppre+imppost) )' * dt;
    
    if ( size( signal_full, 2 ) == 1 )
        signal_impacts = PullImpacts( signal_full, event_inds, imppre, imppost );
        plot( t, signal_impacts );
        plot( t, mean( signal_impacts, 2 ), 'k', 'LineWidth', 2 );
        title( ['Event Traces (n = ', int2str( length( event_inds ) ), ')' ] );
        
        func_out.signal_impacts = signal_impacts;
    else
        assert( size( signal_full, 2 ) == 3 );
        signal_x_impacts = PullImpacts( signal_full(:,1), event_inds, imppre, imppost );
        signal_y_impacts = PullImpacts( signal_full(:,2), event_inds, imppre, imppost );
        signal_z_impacts = PullImpacts( signal_full(:,3), event_inds, imppre, imppost );
        signal_mag_impacts = PullImpacts( sqrt( sum( signal_full .* signal_full, 2 ) ), event_inds, imppre, imppost );
        
        subplot(2,2,1); hold on;
        plot( t, signal_x_impacts );
        plot( t, mean( signal_x_impacts, 2 ), 'k', 'LineWidth', 2 );
        title( ['Event Traces X (n = ', int2str( length( event_inds ) ), ')' ] );
        xlabel( 'Time (s)' )
        set(gca,'FontSize',20)
        
        subplot(2,2,2); hold on;
        plot( t, signal_y_impacts );
        plot( t, mean( signal_y_impacts, 2 ), 'k', 'LineWidth', 2 );
        title( ['Event Traces Y (n = ', int2str( length( event_inds ) ), ')' ] );
        xlabel( 'Time (s)' )
        set(gca,'FontSize',20)
        
        subplot(2,2,3); hold on;
        plot( t, signal_z_impacts );
        plot( t, mean( signal_z_impacts, 2 ), 'k', 'LineWidth', 2 );
        title( ['Event Traces Z (n = ', int2str( length( event_inds ) ), ')' ] );
        xlabel( 'Time (s)' )
        set(gca,'FontSize',20)
        
        subplot(2,2,4); hold on;
        plot( t, signal_mag_impacts );
        plot( t, mean( signal_mag_impacts, 2 ), 'k', 'LineWidth', 2 );
        title( ['Event Traces Mag (n = ', int2str( length( event_inds ) ), ')' ] );
        xlabel( 'Time (s)' )
        set(gca,'FontSize',20)
        
        func_out.signal_x_impacts = signal_x_impacts;
        func_out.signal_y_impacts = signal_y_impacts;
        func_out.signal_z_impacts = signal_z_impacts;
        func_out.signal_mag_impacts = signal_mag_impacts;
    end
    
end