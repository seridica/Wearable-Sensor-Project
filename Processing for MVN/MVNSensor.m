%%%
% File name: MVNSensor.m
% Author: Calvin Kuo
% Date: 11/14/2018
% 
% Extract number of sensors and all sensor data into a n*(m+1) matrix (n 
% time stamps, m sensors and a time column). Updated for MVN 2019, which no
% longer stores the raw sensor data (stores orientation as 4-quaternion and
% 3-axis free acceleration data)

function [sensorNames, sensorData] = MVNSensor(mvnStruct)
    % Get number of sensors
    numSensors = length( mvnStruct.mvnx.subject.sensors.sensor );

    % Get sensor names and fill
    sensorNames = {};
    for i=1:numSensors
        
        % Get current sensor name
        currName = mvnStruct.mvnx.subject.sensors.sensor{i}.Attributes;
        sensorNames = [sensorNames, currName];
    end
    
    % Get number of frames
    nFrames = length( mvnStruct.mvnx.subject.frames.frame );
    
    % Fill sensor data matrix
    sensorData = zeros(nFrames-3, (numSensors*7)+1);
    for i=4:nFrames
        % Time
        sensorData(i-3,1) = str2double( mvnStruct.mvnx.subject.frames.frame{i}.Attributes.time ) / 1000.0;
        
        % Get values to extract
        sensorAccelerations = strsplit( mvnStruct.mvnx.subject.frames.frame{i}.sensorFreeAcceleration.Text, ' ');
        sensorOrientation = strsplit( mvnStruct.mvnx.subject.frames.frame{i}.sensorOrientation.Text, ' ');
        
        % Sensors
        for j=1:numSensors
            sensorData(i-3,(j-1)*7+1+1) = str2double( sensorOrientation{(j-1)*4+1} );
            sensorData(i-3,(j-1)*7+2+1) = str2double( sensorOrientation{(j-1)*4+2} );
            sensorData(i-3,(j-1)*7+3+1) = str2double( sensorOrientation{(j-1)*4+3} );
            sensorData(i-3,(j-1)*7+4+1) = str2double( sensorOrientation{(j-1)*4+4} );
            sensorData(i-3,(j-1)*7+5+1) = str2double( sensorAccelerations{(j-1)*3+1} );
            sensorData(i-3,(j-1)*7+6+1) = str2double( sensorAccelerations{(j-1)*3+2} );
            sensorData(i-3,(j-1)*7+7+1) = str2double( sensorAccelerations{(j-1)*3+3} );
        end
    end
end