%% Helper function to find peaks in a subsection of the full signal
function peak_inds = FindPeaks( t_full, signal_full, t_start, t_end, signal_thresh, wrange )
    
    % Pull section of signal and time
    tinit = find( t_full > t_start, 1, 'first' );
    tfin = find( t_full > t_end, 1, 'first' );
    
    % Make magnitude if not already magnitude
    if size( signal_full, 2 ) == 1
        signal_section = signal_full( tinit:tfin );
    else
        signal_section = sqrt( sum( signal_full( tinit:tfin, : ) .* signal_full( tinit:tfin, : ), 2 ) );
    end
    
    % Cycle through signal to find peaks
    curr_start_ind = 1;
    peak_inds = [];
    peak_ind = find( signal_section( curr_start_ind:end ) > signal_thresh, 1, 'first' );
    while( ~isempty( peak_ind ) )
        
        % Find impact peak
        window_init = curr_start_ind + peak_ind - wrange;
        window_fin = curr_start_ind + peak_ind + wrange;
        %mind = find( signal_section( window_init:window_fin ) > 8, 1, 'first' );
        mind = find( signal_section( window_init:window_fin ) > signal_thresh, 1, 'first' );
        %[p, mind] = max( signal_section( window_init:window_fin) );
        peak_inds = [peak_inds, mind + tinit + window_init];
        %peak_inds = [peak_inds, tinit + curr_start_ind + peak_ind - 1];
        
        % Look for next peak
        curr_start_ind = curr_start_ind + mind + wrange;
        %curr_start_ind = curr_start_ind + peak_ind + wrange;
        peak_ind = find( signal_section( curr_start_ind:end ) > signal_thresh, 1, 'first' );
    end

end