%% Header
% AA228 Final Project
% Calvin Kuo
% 12/6/2016
% UnscentedFilter.m

%% File Description
% This file runs a Unscented filter on three sensors to obtain a
% better estimate of the head impact kinematics. Since the dynamics are
% linear and each sensor is assumed to have Gaussian noise, we can take
% advantage of most of the Kalman filter framework. The only nonlinear
% aspect of our problem getting angular acceleration from a linear
% acceleration network. For this nonlinear portion, we will use the
% unscented framework to generate particle samples and then reconstruct the
% gaussian estimate of the angular acceleration (basically estimate angular
% acceleration and the process noise).

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

function filtered_data = UnscentedFilter( patch_1, patch_2, patch_3, patch_1T, patch_2T, patch_3T, patch_1std, patch_2std, patch_3std )
    
    % Flag for debugging
    debug = 0;

    % Pull useful values out
    patch1_linacc_std = eye(3);
    patch1_linacc_std(1,1) = patch_1std(1,1);
    patch1_linacc_std(2,2) = patch_1std(1,2);
    patch1_linacc_std(3,3) = patch_1std(1,3);
    
    patch2_linacc_std = eye(3);
    patch2_linacc_std(1,1) = patch_2std(1,1);
    patch2_linacc_std(2,2) = patch_2std(1,2);
    patch2_linacc_std(3,3) = patch_2std(1,3);
    
    patch3_linacc_std = eye(3);
    patch3_linacc_std(1,1) = patch_3std(1,1);
    patch3_linacc_std(2,2) = patch_3std(1,2);
    patch3_linacc_std(3,3) = patch_3std(1,3);
    
    patch1_angvel_std = eye(3);
    patch1_angvel_std(1,1) = patch_1std(2,1);
    patch1_angvel_std(2,2) = patch_1std(2,2);
    patch1_angvel_std(3,3) = patch_1std(2,3);
    
    patch2_angvel_std = eye(3);
    patch2_angvel_std(1,1) = patch_2std(2,1);
    patch2_angvel_std(2,2) = patch_2std(2,2);
    patch2_angvel_std(3,3) = patch_2std(2,3);
    
    patch3_angvel_std = eye(3);
    patch3_angvel_std(1,1) = patch_3std(2,1);
    patch3_angvel_std(2,2) = patch_3std(2,2);
    patch3_angvel_std(3,3) = patch_3std(2,3);
    
    % Ensure all sensors have the same length time vector
    t_vec = patch_1.t;
    assert( length( t_vec ) == length( patch_2.t ) );
    assert( length( t_vec ) == length( patch_3.t ) );
    
    % State Vector
    ang_vel_filt = zeros(size(patch_1.lin_acc));
    ang_vel_std = zeros(3, 3, size(patch_1.lin_acc,1));
    ang_acc_filt = zeros(size(patch_1.lin_acc));
    
    % Initialize Particles
    n_particles = 100;
    angacc_particles = zeros(n_particles,3);
    
    % Initial Estimates
    init_angvel_p1 = GetCurrAngVel( patch_1, patch_1T, 1 );
    init_angvel_p2 = GetCurrAngVel( patch_2, patch_2T, 1 );
    init_angvel_p3 = GetCurrAngVel( patch_3, patch_3T, 1 );
    
    % Fuse initial measurements to get initial estimate and covariance
    Hk = [eye(3); eye(3); eye(3)];
    Rsens = eye(9);
    Rsens(1:3,1:3) = patch1_angvel_std;
    Rsens(4:6,4:6) = patch2_angvel_std;
    Rsens(7:9,7:9) = patch3_angvel_std;
    state_est = mean( [init_angvel_p1', init_angvel_p2', init_angvel_p3'], 2 );
    ang_vel_filt(1,:) = state_est';
    cov_est = diag( ( init_angvel_p1' - state_est ).^2 + ...
        ( init_angvel_p2' - state_est ).^2 + ...
        ( init_angvel_p3' - state_est ).^2 ) / 3;
    ang_vel_std(:,:,1) = cov_est;
    
    %% Kalman Filter
    for i=2:length( t_vec )
        % Propagate State and Covariance
        dt = t_vec(i) - t_vec(i-1);
        
        lin_1 = patch_1.lin_acc(i-1,:);
        lin_2 = patch_2.lin_acc(i-1,:);
        lin_3 = patch_3.lin_acc(i-1,:);
            
        % Sample angular acceleation space
        ang_acc_particles = zeros( n_particles, 3 );
        w_vec = zeros( n_particles, 1 );
        for j=1:n_particles
            lin_1_sample = mvnrnd( lin_1, patch1_linacc_std );
            lin_2_sample = mvnrnd( lin_2, patch2_linacc_std );
            lin_3_sample = mvnrnd( lin_3, patch3_linacc_std );
            ang_vel_sample = mvnrnd( state_est, cov_est );
            
            % Weighted based on the probability distribution
            wi = 1;
            wi = wi*mvnpdf(lin_1_sample,lin_1,patch1_linacc_std);
            wi = wi*mvnpdf(lin_2_sample,lin_2,patch2_linacc_std);
            wi = wi*mvnpdf(lin_3_sample,lin_3,patch3_linacc_std);
            w_vec(j) = wi*mvnpdf(ang_vel_sample,state_est',cov_est);
            
            switch randi([1,3])
                case 1
                    lin_a = lin_1_sample'; lin_aT = patch_1T;
                    lin_b = lin_2_sample'; lin_bT = patch_2T;
                    lin_c = lin_3_sample'; lin_cT = patch_3T;
                case 2
                    lin_c = lin_1_sample'; lin_cT = patch_1T;
                    lin_a = lin_2_sample'; lin_aT = patch_2T;
                    lin_b = lin_3_sample'; lin_bT = patch_3T;
                case 3
                    lin_b = lin_1_sample'; lin_bT = patch_1T;
                    lin_c = lin_2_sample'; lin_cT = patch_2T;
                    lin_a = lin_3_sample'; lin_aT = patch_3T;
            end
            
            ang_acc_particles(j,:) = ComputeAngAccel( ang_vel_sample', lin_a, lin_b, lin_c, lin_aT, lin_bT, lin_cT );
        end
        w_vec = w_vec / sum( w_vec );
        ang_acc = w_vec' * ang_acc_particles;
        
        Qsim = zeros(3,3);
        for j=1:3
            %for k=1:3
                %Qsim(j,k) = sum( w_vec.*( ang_acc_particles(:,j) - ang_acc(j) ) .* ( ang_acc_particles(:,k) - ang_acc(k) ) );
            %end
            Qsim(j,j) = sum( w_vec.*( ang_acc_particles(:,j) - ang_acc(j) ) .* ( ang_acc_particles(:,j) - ang_acc(j) ) );
        end
        
        state_est = ang_acc' * dt + state_est;
        cov_est = cov_est + Qsim * dt^2;
        
        % Observations
        angvel_p1 = GetCurrAngVel( patch_1, patch_1T, i );
        angvel_p2 = GetCurrAngVel( patch_2, patch_2T, i );
        angvel_p3 = GetCurrAngVel( patch_3, patch_3T, i );
        
        % Fuse all estimates of the new state and covariance
        Hk = [ eye(3); eye(3); eye(3) ];
        yest = [ angvel_p1'; angvel_p2'; angvel_p3' ] - Hk*state_est;
        S = Hk*cov_est*Hk' + Rsens;
        Kgain = cov_est*Hk'/S;
        state_est = state_est + Kgain * yest;
        cov_est = ( eye(3) - Kgain*Hk ) * cov_est;
        ang_vel_filt(i,:) = state_est';
        ang_vel_std(:,:,i) = cov_est;
        
        if debug == 1
            figure(1); clf;
            curr_t = t_vec(i);
            subplot(1,3,1); hold on;
            plot( ones(n_particles,1)*curr_t, angacc_particles(:,1), 'bx' )
            plot( [curr_t, curr_t], [state_est(1) - cov_est(1,1), state_est(1) + cov_est(1,1)], 'k-' );
            plot( t_vec(1:i), ang_vel_filt(1:i,1), 'r+' )
            xlim([t_vec(1), t_vec(end)])

            subplot(1,3,2); hold on;
            plot( ones(n_particles,1)*curr_t, angacc_particles(:,2), 'bx' )
            plot( [curr_t, curr_t], [state_est(2) - cov_est(2,2), state_est(2) + cov_est(2,2)], 'k-' );
            plot( t_vec(1:i), ang_vel_filt(1:i,2), 'r+' )
            xlim([t_vec(1), t_vec(end)])

            subplot(1,3,3); hold on;
            plot( ones(n_particles,1)*curr_t, angacc_particles(:,3), 'bx' )
            plot( [curr_t, curr_t], [state_est(3) - cov_est(3,3), state_est(3) + cov_est(3,3)], 'k-' );
            plot( t_vec(1:i), ang_vel_filt(1:i,3), 'r+' )
            xlim([t_vec(1), t_vec(end)])
            keyboard;
        end
    end
    
    %% Based on velocity, recalculate other state variables
    % Angular acceleration is derivative of the angular velocity
    
    % Linear acceleration is a weighted average of all the measured linear
    % accelerations.
    
    filtered_data.t = t_vec;
    filtered_data.ang_vel = ang_vel_filt;
    filtered_data.ang_acc = ang_acc_filt;
    filtered_data.ang_vel_std = ang_vel_std;
    %filtered_data.lin_acc = lin_acc_filt;
    %filtered_data.ang_acc = ang_acc_filt;
end

%% Helper function that obtains the angular velocity
% Angular velocity at t_curr and rotated into the head frame.
function curr_angvel = GetCurrAngVel( indata, inT, curr_t )
    this_angvel = indata.ang_vel(curr_t,:);
    R = inT(1:3,1:3);
    curr_angvel = this_angvel*R;
end