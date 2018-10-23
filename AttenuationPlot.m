%% Helper function for attenuation
function AttenuationPlot( signal_orig, Fs, freqs, unique_inds, imppre, imppost )

    % Set up the peak and attenuation matrix
    attenuation_matrix = zeros( length( freqs )+1, length( unique_inds ) );
    peak_matrix = zeros( length( freqs )+1, length( unique_inds ) );
    
    % Find peaks for the input signal
    for i=1:length( unique_inds )
        peak_ind = unique_inds(i);
        impact_signal = signal_orig(peak_ind-imppre:peak_ind+imppost);
        peak_matrix(1,i) = max( abs( impact_signal ) );
    end
    
    % Go through filters
    for i=1:length( freqs )
        [b,a] = butter( 2, freqs(i) / (Fs/2) );
        filt_signal = filtfilt( b, a, signal_orig );
        for j=1:length( unique_inds )
            peak_ind = unique_inds(j);
            impact_signal = filt_signal(peak_ind-imppre:peak_ind+imppost);
            peak_matrix(i+1,j) = max( abs( impact_signal ) );
        end
        attenuation_matrix(i+1,:) = ( peak_matrix(1,:) - peak_matrix(i+1,:) ) ./ peak_matrix(1,:);
    end
    
    subplot(2,1,1);
    notBoxPlot( attenuation_matrix' );
    ylabel( 'Attenuation' )
    xlabel( 'Filter (Hz' )
    set( gca, 'FontSize', 20 );
    
    subplot(2,1,2);
    notBoxPlot( peak_matrix' );
    ylabel( 'Peak' )
    xlabel( 'Filter (Hz' )
    set( gca, 'FontSize', 20 );
end