%%%
% File name: MVNJointAngle.m
% Author: Calvin Kuo
% Date: 11/14/2018
%
% From an mvnstruct file, pull the joint angles
% information

function [jointNames, jointData] = MVNJointAngle(mvnStruct)
    % Get number of sensors
    numJoints = length( mvnStruct.mvnx.subject.joints.joint );

    % Get sensor names and fill
    jointNames = {};
    for i=1:numJoints
        
        % Get current sensor name
        currName = mvnStruct.mvnx.subject.segments.segment{i}.Attributes;
        jointNames = [jointNames, currName];
    end
    
    % Get number of frames
    nFrames = length( mvnStruct.mvnx.subject.frames.frame );
    
    % Fill sensor data matrix
    jointData = zeros(nFrames-3, (numJoints*3)+1);
    for i=4:nFrames
        % Time
        jointData(i-3,1) = str2double( mvnStruct.mvnx.subject.frames.frame{i}.Attributes.time ) / 1000.0;
        
        % Get values to extract
        jointAngles = strsplit( mvnStruct.mvnx.subject.frames.frame{i}.jointAngle.Text, ' ');
        
        % Sensors
        for j=1:numJoints
            jointData(i-3,(j-1)*3+1+1) = str2double( jointAngles{(j-1)*3+1} );
            jointData(i-3,(j-1)*3+2+1) = str2double( jointAngles{(j-1)*3+2} );
            jointData(i-3,(j-1)*3+3+1) = str2double( jointAngles{(j-1)*3+3} );
        end
    end
end