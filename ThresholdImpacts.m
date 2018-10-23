%%%
% File name: FindImpacts.m
% Author: Calvin Kuo
% Date: 9/21/2018
% 
% Parse linear acceleratio traces to look for impacts in the signal.
% Returns a structure of the impacts containing the time in which they
% occurred, index range from the source, and the data (time + linear
% accelerations). Uses standard thresholding to look for impacts (over
% 20m/s^2 in the magnitude linear acceleration, -100ms, +800ms)

function impacts = ThresholdImpacts(tvec, linacc, threshold, pre, post)

    % Set up impacts list
    impacts = [];
    
    % Parse the time vector for dt and length of signal
    nSamples = length (tvec);
    dt = ( tvec(end) - tvec(1) ) / nSamples;
    preInds = floor( pre / dt );
    postInds = floor( post / dt );
    
    % Magnitude linear acceleration
    magLin = sqrt( sum( linacc .* linacc, 2 ) );
    
    % Loop through signal to look for impacts
    i = preInds;
    while ( i < (nSamples-postInds) )
    %for i=preInds:(nSamples-postInds)
        
        % Check threshold
        if (magLin(i) >= threshold)
            
            % Make an impact structure
            impStruct = struct();
            impStruct.impact_time = tvec(i);
            impStruct.t = tvec((i-preInds):(i+postInds)) - tvec(i);
            impStruct.linaccMag = magLin((i-preInds):(i+postInds));
            impStruct.linacc = linacc((i-preInds):(i+postInds),:);
            
            impacts = [impacts, impStruct];
            
            % Move time signal forward past impact
            i = i+postInds;
        end
        i = i+1;
    end
    
end