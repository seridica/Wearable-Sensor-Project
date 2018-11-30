%% Header
% AA228 Final Project
% Calvin Kuo
% 12/6/2016
% KalmanFilterNew.m

%% File Description
% This file runs a Kalman filter on three sensors to obtain a
% better estimate of the head impact kinematics. Since the dynamics are
% linear and each sensor is assumed to have Gaussian noise, we can take
% advantage of most of the Kalman filter framework. The only nonlinear
% aspect of our problem getting angular acceleration from a linear
% acceleration network, but we assume the nonlinear term (omega cross
% omega), can be ignored.

%% Function
% Inputs
% 1) Patch data in the form of structs containing UNTRANSFORMED triaxial
% linear acceleration and angular velocity measurements and the time vector
% 2) Patch transforms for each patch to a reference.
% 3) Estimated standard deviation for the patch sensor measurement errors
%
% Outputs
% 1) Filtered data struct containing the time vector, the estimated
% triaxial linear acceleration and angular velocity, as well as the
% estimated standard deviations for the estimtes.

function filtered_data = KalmanFilterNew( patch_data, patch_transforms, patch_stds, std_scale )

    % Scaling error
    if nargin == 3
        std_scale = 0;
    end
    
    % Flag for debugging
    debug = 0;
        
    assert( length( patch_data ) == length( patch_transforms ) );
    assert( length( patch_data ) == length( patch_stds ) );
    
    ndat = length( patch_data );
    
    % Pull useful values out
    patch_linacc_std = zeros(ndat,3);
    patch_angvel_std = zeros(ndat,3);
    for i=1:ndat
        patch_linacc_std(i,:) = patch_stds{i}(1,:);
        patch_angvel_std(i,:) = patch_stds{i}(2,:);
    end
    
    % Ensure all sensors have the same length time vector
    t_vec = patch_data{1}.t;
    for i=1:ndat
        t_check = patch_data{i}.t;
        assert( length( t_vec ) == length( t_check ) );
    end
    
    % Average patch angular velocities
    filtered_data.ang_vel = zeros( length( t_vec ), 3 );
    
    % Data Structures
    filtered_data.ang_acc = zeros( length( t_vec ), 3 );
    filtered_data.lin_acc = zeros( length( t_vec ), 3 );
    
    %% Particle filter to determine states
    R_v_base = zeros(ndat*6,ndat*6);
    for i=1:ndat
        R_v_base(((i-1)*6+1):((i-1)*6+3), ((i-1)*6+1):((i-1)*6+3)) = diag( patch_angvel_std(i,:) );
        R_v_base(((i-1)*6+4):(i*6), ((i-1)*6+4):(i*6)) = diag( patch_linacc_std(i,:) );
    end
    
    if ( std_scale == 1 )
        R_v = R_v_base;
    else
        R_v = R_v_base;
    end
    
    % Build H matrix
    H = [];
    for i=1:ndat
        this_t = patch_transforms{i}(1:3,4) * 2.54 / 100;
        cross_mat = [0, this_t(3), -this_t(2); ...
            -this_t(3), 0, this_t(1); ...
            this_t(2), -this_t(1), 0];
        this_H = [eye(3), zeros(3,3), zeros(3,3); ...
            zeros(3,3), cross_mat, eye(3)];
        H = [H; this_H];
    end
    
    % First estimate
    mvec = [];
    for i=1:ndat
        this_m = ( GetCurrMeas( patch_data{i}, patch_transforms{i}, 1 ) )';
        mvec = [mvec; this_m];
    end
    
    curr_cov = inv(H'*inv(R_v)*H);
    curr_state = curr_cov*H'*inv(R_v)*mvec;
    %keyboard;

    filtered_data.ang_vel(1,:) = curr_state(1:3)';
    filtered_data.ang_acc(1,:) = curr_state(4:6)';
    filtered_data.lin_acc(1,:) = curr_state(7:9)' / 9.81;

    if ( std_scale == 1 )
        mav = norm( filtered_data.ang_vel(1,:), 2 );
        mla = norm( filtered_data.lin_acc(1,:), 2 );

        full_temp = zeros( 6*ndat, 6*ndat );

        for j=1:ndat
            dtemp = [max( patch_angvel_std(j,:)*mav, 0.1), max( patch_linacc_std(j,:)*mla, 0.5 )];
            temp = diag( dtemp );
            full_temp( ((j-1)*6+1):(j*6), ((j-1)*6+1):(j*6) ) = temp;
        end
        
        R_v = R_v_base * full_temp;
    else
        R_v = R_v_base;
    end
    
    for i=2:length( t_vec )
        % State Propagation
        dt = t_vec(i) - t_vec(i-1);
        F = zeros(9,9);
        F(1:3,1:3) = eye(3);
        F(1:3,4:6) = eye(3) * dt;
        
        apri_state = F*curr_state;
        apri_cov = F*curr_cov*F';
        
        inv_cov = zeros(9,9);
        inv_cov(1:3, 1:3) = inv(apri_cov(1:3,1:3));
        
        mvec = [];
        for j=1:ndat
            this_m = ( GetCurrMeas( patch_data{j}, patch_transforms{j}, i ) )';
            mvec = [mvec; this_m];
        end

        yest = mvec - H*apri_state;
        curr_cov = inv(inv_cov + H'*inv(R_v)*H);
        curr_state = apri_state + curr_cov*H'*inv(R_v)*yest;
        
        %S = H*apri_cov*H' + R_v;
        %Kgain = apri_cov*H'/S;
        %curr_state = apri_state + Kgain*[m1'; m2'; m3'];
        %curr_cov = ( eye(9) - Kgain*H ) * apri_cov;
    
        %keyboard;
        
        filtered_data.ang_vel(i,:) = curr_state(1:3)';
        filtered_data.ang_acc(i,:) = curr_state(4:6)';
        filtered_data.lin_acc(i,:) = curr_state(7:9)' / 9.81;
        
        if ( std_scale == 1 )
            mav = norm( filtered_data.ang_vel(i,:), 2 );
            mla = norm( filtered_data.lin_acc(i,:), 2 );

            full_temp = zeros( 6*ndat, 6*ndat );

            for j=1:ndat
                dtemp = [max( patch_angvel_std(j,:)*mav, 0.1), max( patch_linacc_std(j,:)*mla, 0.5 )];
                temp = diag( dtemp );
                full_temp( ((j-1)*6+1):(j*6), ((j-1)*6+1):(j*6) ) = temp;
            end
        
            R_v = R_v_base * full_temp;
        else
            R_v = R_v_base;
        end
    end
    
    %% Based on velocity, recalculate other state variables
    % Angular acceleration is derivative of the angular velocity
    
    % Linear acceleration is a weighted average of all the measured linear
    % accelerations.
    
    filtered_data.t = t_vec;
    filtered_data.ang_vel_std = zeros( length( t_vec ), 3 );
    filtered_data.ang_acc_std = zeros( length( t_vec ), 3 );
    filtered_data.lin_acc_std = zeros( length( t_vec ), 3 );
end

%% Helper function that obtains the current measurement vector
% measurement vector at time t
function curr_meas = GetCurrMeas( indata, inT, curr_t )
    this_angvel = indata.ang_vel(curr_t,:);
    this_linacc = indata.lin_acc(curr_t,:);
    R = inT(1:3,1:3);
    curr_angvel = this_angvel*R;
    curr_linacc = this_linacc*R*9.81;
    
    curr_meas = [curr_angvel, curr_linacc];
end