%% Helper function to plot fft
function HelperPlotFft( t_full, signal_full, t_start, t_end, freq_cutoff )
    
    % Pull section of signal and time
    tinit = find( t_full > t_start, 1, 'first' );
    tfin = find( t_full > t_end, 1, 'first' );
    
    t = t_full( tinit:tfin );        % Time vector
    L = length( t );
    T = ( t(end) - t(1) ) / L;
    Fs = 1/T;
    
    signal_section = signal_full( tinit:tfin );
    
    % Fft from matlab documentation - Demean the signal
    Y = fft(signal_section - mean( signal_section ));
    
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    subplot( 2, 1, 1 );
    plot( t, signal_section );
    
    subplot( 2, 1, 2 );
    plot(f,P1)
    xlim([0, freq_cutoff]);
    
end