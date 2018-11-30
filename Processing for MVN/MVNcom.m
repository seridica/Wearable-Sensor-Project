%%%
% File name: MVNcom.m
% Author: Calvin Kuo
% Date: 11/14/2018
%
% From an mvnstruct file, pull the center of mass information

function comPos = MVNcom(mvnStruct)
    
    % Get number of frames
    nFrames = length( mvnStruct.mvnx.subject.frames.frame );
    
    % Fill sensor data matrix
    comPos = zeros(nFrames-3, 4);
    for i=4:nFrames
        % Time
        comPos(i-3,1) = str2double( mvnStruct.mvnx.subject.frames.frame{i}.Attributes.time ) / 1000.0;
        
        % Get values to extract
        comPosition = strsplit( mvnStruct.mvnx.subject.frames.frame{i}.centerOfMass.Text, ' ');
        
        % Sensors
        comPos(i-3,2:4) = str2double( comPosition );
    end
end