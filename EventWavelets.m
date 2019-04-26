%% Helper function for wavelet transforms of the impacts
function maxabs = EventWavelets( Fs, event_signals, peak_percent, crange )

    % Cycle through event signals for the wavelet transform
    nEvents = size( event_signals, 2 );
    signal_cfs = [];
    for i=1:nEvents
        this_signal = event_signals(:,i);
        [cfs, f] = cwt( this_signal, Fs );
        nRows = size( cfs, 2 );
        nCols = size( cfs, 1 );
        signal_cfs = [signal_cfs, abs( reshape( cfs, nRows*nCols, 1 ) )];
    end
    
    t = ( 0:( size( event_signals, 1 ) - 1 ) ) / Fs;
    surface( 'XData',t,'YData',f,...
         'CData',reshape( mean(signal_cfs')', nCols, nRows), 'ZData', zeros(nCols,nRows), ...
         'CDataMapping','scaled', 'FaceColor','texturemap', 'EdgeColor', 'none');
    ylim([0,5]);
    colorbar();
    if (nargin > 3)
        caxis(crange)
    end
    
    % Contour overlay
    % 5%, 10%, 20%
    [M, c] = contour( t, f, reshape( mean(signal_cfs')', nCols, nRows), [0.2, 0.2]*peak_percent, 'k:' );
    c.LineWidth = 2.0;
    [M, c] = contour( t, f, reshape( mean(signal_cfs')', nCols, nRows), [0.1, 0.1]*peak_percent, 'k-.' );
    c.LineWidth = 2.0;
    [M, c] = contour( t, f, reshape( mean(signal_cfs')', nCols, nRows), [0.05, 0.05]*peak_percent, 'k' );
    c.LineWidth = 2.0;
    
    maxabs = max( mean( signal_cfs' )' );
end