%%%
% File name: FindTransform.m
% Author: Calvin Kuo
% Date: 10/19/2018
%
% Finds the sensor transform given three points on a sensor
% Point 1 - sensor frame origin
% Point 2 - along x-axis
% Point 3 - along y-axis
%
% Returns transform containing rotation from lab frame to sensor frame 
% defined as above, and the translation from lab frame origin to point 1.

function transform_matrix = FindTransform( p1, p2, p3 )
    x_axis_long = p2 - p1;
    x_axis = x_axis_long / norm( x_axis_long );
    
    y_axis_temp_long = p3 - p1;
    y_axis_temp = y_axis_temp_long / norm( y_axis_temp_long );
    
    z_axis_long = cross( x_axis, y_axis_temp );
    z_axis = z_axis_long / norm( z_axis_long );
    
    y_axis_long = cross( z_axis, x_axis );
    y_axis = y_axis_long / norm( y_axis_long );
    
    transform_matrix = eye(4);
    transform_matrix(1:3,4) = p1;
    transform_matrix(1:3,1) = x_axis;
    transform_matrix(1:3,2) = y_axis;
    transform_matrix(1:3,3) = z_axis;
%     transform_matrix(1,1:3) = x_axis';
%     transform_matrix(2,1:3) = y_axis';
%     transform_matrix(3,1:3) = z_axis';
end