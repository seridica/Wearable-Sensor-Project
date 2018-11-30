%% Header
% AA228 Final Project
% Calvin Kuo
% 5/12/2016
% ProjectLinAcc.m

%% File Description
% Projecs linear acceleration given lin acc (in sensor frame)
% ang acc, ang vel in world frame, and transformation from world to sensor.

%% Function
function proj_acc = ProjectLinAcc( lin_acc, ang_vel, ang_acc, T )
    g = 9.81;
    r_vec = T(1:3,4) * 2.54 / 100;
    R = T(1:3,1:3);
    wrld_acc = lin_acc'*R*g;
    proj_acc = ( wrld_acc' + cross( ang_vel', cross( ang_vel', -r_vec ) ) + cross( ang_acc', -r_vec ) )' / g;
end