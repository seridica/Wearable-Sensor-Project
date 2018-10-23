%% Helper function for pulling impact time traces
function impact_traces = PullImpacts( signal_full, event_inds, imppre, imppost )

    % Set up output matrix for data
    num_events = length( event_inds );
    impact_traces = zeros( ( imppost + imppre + 1 ), num_events );
    
    % If signal_full is not a vector, then take the magnitude
    if size( signal_full, 2 ) == 1
        signal_pull = signal_full;
    else
        signal_pull = sqrt( sum( signal_full .* signal_full, 2 ) );
    end
    
    % Pull traces
    for i=1:num_events
        raw_trace = signal_pull( ( event_inds(i)-imppre ):( event_inds(i)+imppost ) );
        %impact_traces(:, i) = raw_trace / ( max( abs( raw_trace ) ) ); % Normalize
        impact_traces(:,i) = raw_trace;
    end
end