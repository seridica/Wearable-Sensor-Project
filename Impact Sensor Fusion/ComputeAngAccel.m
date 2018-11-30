%% Header
% AA228 Final Project
% Calvin Kuo
% 11/21/2016
% ComputeAngAccel.m

%% File Description
% Function that computes angular acceleration calculated from three
% triaxial linear accelerometer blocks arranged arbitrarily (not colinear).
% For the algorithm, see project.

%% Function
% Inputs:
% 1 - Ang Vel Estimate (3-vector)
% 2 - Lin Acc block 1 Estimate (3-vector)
% 3 - Lin Acc block 2 Estimate (3-vector)
% 4 - Lin Acc block 3 Estimate (3-vector)
% 5 - Lin Acc block 1 Transform (4x4 matrix)
% 6 - Lin Acc block 2 Transform (4x4 matrix)
% 7 - Lin Acc block 3 Transform (4x4 matrix)

% Outputs:
% 1 - Ang Acc Estimate (3-axis)

function ang_acc = ComputeAngAccel( ang_vel, lin_acc_1, lin_acc_2, lin_acc_3, T_1, T_2, T_3 )
    % Extract useful things
    g = 9.81;
    p1 = T_1(1:3,4) * 2.54 / 100;
    p2 = T_2(1:3,4) * 2.54 / 100;
    p3 = T_3(1:3,4) * 2.54 / 100;
    
    R1 = T_1(1:3,1:3)';% * [1 0 0; 0 -1 0; 0 0 -1];
    R2 = T_2(1:3,1:3)';% * [1 0 0; 0 -1 0; 0 0 -1];
    R3 = T_3(1:3,1:3)';% * [1 0 0; 0 -1 0; 0 0 -1];
    
    % Find common frame based on position of triaxials
    x_vec = p2-p1;
    x_vec = x_vec / norm( x_vec, 2 );
    
    y_tmp = p3-p1;
    y_tmp = y_tmp / norm( y_tmp, 2 );
    
    z_vec = cross( x_vec, y_tmp );
    z_vec = z_vec / norm( z_vec, 2 );
    
    y_vec = cross( z_vec, x_vec );
    y_vec = y_vec / norm( y_vec, 2 );
    
    R_common = [x_vec, y_vec, z_vec];
    p_common = dot( p3-p1, x_vec ) * x_vec + p1;
    
    % Rotate data to frame
    la1_common = R_common' * R1 * lin_acc_1 * g;
    la2_common = R_common' * R2 * lin_acc_2 * g;
    la3_common = R_common' * R3 * lin_acc_3 * g;
    
    av_common = R_common' * ang_vel;
    
    % Compute angular acceleration
    rx = dot( ( p2 - p1 ), x_vec );
    aa_y = -(la2_common(3) - la1_common(3) - av_common(1)*av_common(3)*rx) / rx;
    aa_z = (la2_common(2) - la1_common(2) - av_common(1)*av_common(2)*rx) / rx;
    
    rx_mid = dot( ( p_common - p1 ), x_vec );
    z_comm = la1_common(3) + av_common(1)*av_common(3)*rx_mid - aa_y*rx_mid;
    ry = dot( ( p3 - p_common ), y_vec );
    aa_x = (la3_common(3) - z_comm - av_common(2)*av_common(3)*ry) / ry;
    
    %aa_common = -[aa_x; aa_y; aa_z];
    %ang_acc = R_common * aa_common;
    %ang_acc = ang_acc([2,1,3]);
    
    aa_common = [aa_x; aa_y; aa_z];
    ang_acc = R_common * aa_common;
end