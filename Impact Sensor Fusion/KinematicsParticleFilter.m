%% Header
% AA228 Final Project
% Calvin Kuo
% 11/26/2016
% KinematicsParticleFilter.m

%% File Description
% This file runs a particle filter on three sensors to obtain a better
% estimate of the head impact kinematics.

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

function filtered_data = KinematicsParticleFilter( patch_1, patch_2, patch_3, patch_1T, patch_2T, patch_3T, patch_1std, patch_2std, patch_3std )
    
    % Flag for debugging
    debug = 0;

    % Pull useful values out
    patch1_linacc_std = patch_1std(1,:);
    patch2_linacc_std = patch_2std(1,:);
    patch3_linacc_std = patch_3std(1,:);
    
    patch1_angvel_std = patch_1std(2,:);
    patch2_angvel_std = patch_2std(2,:);
    patch3_angvel_std = patch_3std(2,:);
    
    % Ensure all sensors have the same length time vector
    t_vec = patch_1.t;
    assert( length( t_vec ) == length( patch_2.t ) );
    assert( length( t_vec ) == length( patch_3.t ) );
    
    % State Vector
    ang_vel_filt = zeros(size(patch_1.lin_acc));
    ang_vel_std = zeros(size(patch_1.lin_acc));
    ang_acc_filt = zeros(size(patch_1.lin_acc));
    lin_acc_filt = zeros(size(patch_1.lin_acc));
    
    % Initial Estimates
    init_angvel_p1 = GetCurrAngVel( patch_1, patch_1T, 1 );
    init_angvel_p2 = GetCurrAngVel( patch_2, patch_2T, 1 );
    init_angvel_p3 = GetCurrAngVel( patch_3, patch_3T, 1 );
    
    % Initialize Particles
    n_particles = 100;
    
    mean_angvels = mean( [init_angvel_p1; init_angvel_p2; init_angvel_p3], 1 );
    std_angvels = std( [init_angvel_p1; init_angvel_p2; init_angvel_p3], 1 );
    particle_states = zeros( n_particles, 3 );
    other_states = zeros( n_particles, 6 );
    for i=1:n_particles
        particle_states(i,:) = mean_angvels - std_angvels + 2*std_angvels.*rand(1,3);
    end
    
    particle_weights = CalculateWeights( particle_states, [init_angvel_p1; init_angvel_p2; init_angvel_p3], [patch1_angvel_std; patch2_angvel_std; patch3_angvel_std] );
    
    % Check if need to resample particles
    n_eff = 1 / ( particle_weights' * particle_weights );
    if n_eff < n_particles / 2
        particle_states = ResampleParticles(particle_states, particle_weights);
        particle_weights = CalculateWeights( particle_states, [init_angvel_p1; init_angvel_p2; init_angvel_p3], [patch1_angvel_std; patch2_angvel_std; patch3_angvel_std] );
    end
    
    ang_vel_filt(1,:) = particle_weights' * particle_states;
    for i=1:size( ang_vel_filt, 2 )
        ang_vel_std(1,i) = particle_weights' * ( particle_states(:,i) - ang_vel_filt(1,i) ).^2;
    end
    
    if debug == 1
        figure(1); clf;
        curr_t = t_vec(1);
        subplot(1,3,1); hold on;
        plot( ones(n_particles,1)*curr_t, particle_states(:,1), 'bx' )
        plot( [curr_t, curr_t], [ang_vel_filt(1,1) - ang_vel_std(1,1), ang_vel_filt(1,1) + ang_vel_std(1,1)], 'k-' );
        plot( curr_t, ang_vel_filt(1,1), 'r+' )
        xlim([t_vec(1), t_vec(end)])
        
        subplot(1,3,2); hold on;
        plot( ones(n_particles,1)*curr_t, particle_states(:,2), 'bx' )
        plot( [curr_t, curr_t], [ang_vel_filt(1,2) - ang_vel_std(1,2), ang_vel_filt(1,2) + ang_vel_std(1,2)], 'k-' );
        plot( curr_t, ang_vel_filt(1,2), 'r+' )
        xlim([t_vec(1), t_vec(end)])
        
        subplot(1,3,3); hold on;
        plot( ones(n_particles,1)*curr_t, particle_states(:,3), 'bx' )
        plot( [curr_t, curr_t], [ang_vel_filt(1,3) - ang_vel_std(1,3), ang_vel_filt(1,3) + ang_vel_std(1,3)], 'k-' );
        plot( curr_t, ang_vel_filt(1,3), 'r+' )
        xlim([t_vec(1), t_vec(end)])
        keyboard;
    end
    
    %% Particle filter to determine angular velocity
    for i=2:length( t_vec )
        
        % Estimate angular acceleration
        dt = t_vec(i) - t_vec(i-1);
        for j=1:n_particles
            lin_1 = patch_1.lin_acc(i-1,:);
            lin_2 = patch_2.lin_acc(i-1,:);
            lin_3 = patch_3.lin_acc(i-1,:);
            
            if ( std_scale == 1 )
                lin_1_sample = normrnd( lin_1, patch1_linacc_std*sqrt(sum(lin_1.*lin_1)) );
                lin_2_sample = normrnd( lin_2, patch2_linacc_std*sqrt(sum(lin_2.*lin_2)) );
                lin_3_sample = normrnd( lin_3, patch3_linacc_std*sqrt(sum(lin_3.*lin_3)) );
            else
                lin_1_sample = normrnd( lin_1, patch1_linacc_std );
                lin_2_sample = normrnd( lin_2, patch2_linacc_std );
                lin_3_sample = normrnd( lin_3, patch3_linacc_std );
            end
            ang_vel_sample = particle_states(j,:);
            
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
            
            ang_acc = ComputeAngAccel( ang_vel_sample', lin_a, lin_b, lin_c, lin_aT, lin_bT, lin_cT );
            other_states(j,1:3) = ang_acc;
            other_states(j,4:6) = ProjectLinAcc( lin_a, ang_vel_sample, other_states(i-1,1:3), lin_aT );
            ang_acc_filt(i-1,:) = particle_weights' * other_states(:,1:3);
            lin_acc_filt(i-1,:) = particle_weights' * other_states(:,4:6);
            
            % Propagate
            particle_states(j,:) = ang_acc' * dt + ang_vel_sample;
        end
        
        % Observations
        angvel_p1 = GetCurrAngVel( patch_1, patch_1T, i );
        angvel_p2 = GetCurrAngVel( patch_2, patch_2T, i );
        angvel_p3 = GetCurrAngVel( patch_3, patch_3T, i );
        particle_weights = CalculateWeights( particle_states, [angvel_p1; angvel_p2; angvel_p3], [patch1_angvel_std; patch2_angvel_std; patch3_angvel_std] );
        
        % Check if need to resample Particles
        n_eff = 1 / ( particle_weights' * particle_weights );
        if n_eff < n_particles / 2
            particle_states = ResampleParticles(particle_states, particle_weights);
            particle_weights = CalculateWeights( particle_states, [angvel_p1; angvel_p2; angvel_p3], [patch1_angvel_std; patch2_angvel_std; patch3_angvel_std] );
        end
        
        % Estimate
        ang_vel_filt(i,:) = particle_weights' * particle_states;
        for j=1:size( ang_vel_filt, 2 )
            ang_vel_std(i,j) = particle_weights' * ( particle_states(:,j) - ang_vel_filt(i,j) ).^2;
        end
        
        if debug == 1
            figure(1); clf;
            curr_t = t_vec(i);
            subplot(1,3,1); hold on;
            plot( ones(n_particles,1)*curr_t, particle_states(:,1), 'bx' )
            plot( [curr_t, curr_t], [ang_vel_filt(i,1) - ang_vel_std(i,1), ang_vel_filt(i,1) + ang_vel_std(i,1)], 'k-' );
            plot( t_vec(1:i), ang_vel_filt(1:i,1), 'r+' )
            xlim([t_vec(1), t_vec(end)])

            subplot(1,3,2); hold on;
            plot( ones(n_particles,1)*curr_t, particle_states(:,2), 'bx' )
            plot( [curr_t, curr_t], [ang_vel_filt(i,2) - ang_vel_std(i,2), ang_vel_filt(i,2) + ang_vel_std(i,2)], 'k-' );
            plot( t_vec(1:i), ang_vel_filt(1:i,2), 'r+' )
            xlim([t_vec(1), t_vec(end)])

            subplot(1,3,3); hold on;
            plot( ones(n_particles,1)*curr_t, particle_states(:,3), 'bx' )
            plot( [curr_t, curr_t], [ang_vel_filt(i,3) - ang_vel_std(i,3), ang_vel_filt(i,3) + ang_vel_std(i,3)], 'k-' );
            plot( t_vec(1:i), ang_vel_filt(1:i,3), 'r+' )
            xlim([t_vec(1), t_vec(end)])
            keyboard;
        end
    end
    
    lin_1 = patch_1.lin_acc(end,:);
    lin_2 = patch_2.lin_acc(end,:);
    lin_3 = patch_3.lin_acc(end,:);

    if ( std_scale == 1 )
        lin_1_sample = normrnd( lin_1, patch1_linacc_std*sqrt(sum(lin_1.*lin_1)) );
        lin_2_sample = normrnd( lin_2, patch2_linacc_std*sqrt(sum(lin_2.*lin_2)) );
        lin_3_sample = normrnd( lin_3, patch3_linacc_std*sqrt(sum(lin_3.*lin_3)) );
    else
        lin_1_sample = normrnd( lin_1, patch1_linacc_std );
        lin_2_sample = normrnd( lin_2, patch2_linacc_std );
        lin_3_sample = normrnd( lin_3, patch3_linacc_std );
    end
    ang_vel_sample = particle_states(j,:);

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

    ang_acc = ComputeAngAccel( ang_vel_sample', lin_a, lin_b, lin_c, lin_aT, lin_bT, lin_cT );
    other_states(end,1:3) = ang_acc;
    other_states(end,4:6) = ProjectLinAcc( lin_a, ang_vel_sample, other_states(end,1:3), lin_aT );
    ang_acc_filt(i-1,:) = particle_weights' * other_states(:,1:3);
    lin_acc_filt(i-1,:) = particle_weights' * other_states(:,4:6);
    
    %% Based on velocity, recalculate other state variables
    % Angular acceleration is derivative of the angular velocity
    
    % Linear acceleration is a weighted average of all the measured linear
    % accelerations.
    
    filtered_data.t = t_vec;
    filtered_data.ang_vel = ang_vel_filt;
    filtered_data.ang_acc = ang_acc_filt;
    filtered_data.ang_vel_std = ang_vel_std;
    filtered_data.lin_acc = lin_acc_filt;
end

%% Helper function that obtains the angular velocity
% Angular velocity at t_curr and rotated into the head frame.
function curr_angvel = GetCurrAngVel( indata, inT, curr_t )
    this_angvel = indata.ang_vel(curr_t,:);
    R = inT(1:3,1:3);
    curr_angvel = this_angvel*R;
end

%% Helper function that resamples particles based on distribution
% Stratified Resampling
% https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python/blob/master/12-Particle-Filters.ipynb
function new_particles = ResampleParticles( curr_particles, curr_weights )
%    new_particles = curr_particles;
    N = length( curr_weights );
    pos1 = rand(1,N);
    pos2 = 0:N-1;
    positions = (pos1+pos2) / N;
    indexes = zeros(N,1);
    cumulative_sum = cumsum( curr_weights );
    i = 1; j = 1;
    while i<=N
        if positions(i) < cumulative_sum(j);
            indexes(i) = j;
            i = i+1;
        else
            j = j+1;
        end
    end
    new_particles = curr_particles(indexes,:);
end

% Systematic Resampling
% https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python/blob/master/12-Particle-Filters.ipynb