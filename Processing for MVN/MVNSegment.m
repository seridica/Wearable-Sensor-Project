%%%
% File name: MVNSegment.m
% Author: Calvin Kuo
% Date: 11/14/2018
%
% From an mvnstruct file, pull the segment position and orientation
% information

function [segmentNames, segmentData] = MVNSegment(mvnStruct)
    % Get number of sensors
    numSegments = length( mvnStruct.mvnx.subject.segments.segment );

    % Get sensor names and fill
    segmentNames = {};
    for i=1:numSegments
        
        % Get current sensor name
        currName = mvnStruct.mvnx.subject.segments.segment{i}.Attributes;
        segmentNames = [segmentNames, currName];
    end
    
    % Get number of frames
    nFrames = length( mvnStruct.mvnx.subject.frames.frame );
    
    % Fill sensor data matrix
    segmentData = zeros(nFrames-3, (numSegments*6)+1);
    for i=4:nFrames
        % Time
        segmentData(i-3,1) = str2double( mvnStruct.mvnx.subject.frames.frame{i}.Attributes.time ) / 1000.0;
        
        % Get values to extract
        segmentPosition = strsplit( mvnStruct.mvnx.subject.frames.frame{i}.position.Text, ' ');
        segmentOrientation = strsplit( mvnStruct.mvnx.subject.frames.frame{i}.orientation.Text, ' ');
        
        % Sensors
        for j=1:numSegments
            segmentData(i-3,(j-1)*6+1+1) = str2double( segmentOrientation{(j-1)*3+1} );
            segmentData(i-3,(j-1)*6+2+1) = str2double( segmentOrientation{(j-1)*3+2} );
            segmentData(i-3,(j-1)*6+3+1) = str2double( segmentOrientation{(j-1)*3+3} );
            segmentData(i-3,(j-1)*6+4+1) = str2double( segmentPosition{(j-1)*3+1} );
            segmentData(i-3,(j-1)*6+5+1) = str2double( segmentPosition{(j-1)*3+2} );
            segmentData(i-3,(j-1)*6+6+1) = str2double( segmentPosition{(j-1)*3+3} );
        end
    end
end