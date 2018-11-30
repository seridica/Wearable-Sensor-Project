%% Header
% AA228 Final Project
% Calvin Kuo
% 5/11/2016
% ImprovedKinematicsParticleFilter.m

%% File Description
% This file runs a particle filter on n sensors to obtain a better
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

function filtered_data = ImprovedKinematicsParticleFilter( patch_data, patch_transforms, patch_stds, std_scale )
    
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
    
    % State Vector - 9vec state: 3vec ang vel, 3vec ang acc, 3vec lin acc
    nstates = 9;
    state_mat = zeros(length(t_vec),nstates);
    state_std = zeros(length(t_vec),nstates);
    
    % Initial Estimates
    all_angvels = zeros(ndat,3);
    for i=1:ndat
        all_angvels(i,:) = GetCurrAngVel( patch_data{i}, patch_transforms{i}, 1 );
    end
    
    % Initialize Particles
    n_particles = 5000;
    
    mean_angvels = mean( all_angvels, 1 );
    std_angvels = std( all_angvels, 1 );
    particle_states = zeros( n_particles, nstates );
    for i=1:n_particles
        particle_states(i,1:3) = mean_angvels - std_angvels + 2*std_angvels.*rand(1,3);
    
        % Compute angular acceleration state estimate
        sens = randperm(ndat);
        ang_vel_sample = particle_states(i,1:3);
        
        lin_1 = patch_data{sens(1)}.lin_acc(1,:);
        lin_2 = patch_data{sens(2)}.lin_acc(1,:);
        lin_3 = patch_data{sens(3)}.lin_acc(1,:);
        
        if ( std_scale == 1 )
            lin_1_sample = normrnd( lin_1, patch_linacc_std(sens(1))*sqrt(sum(lin_1.*lin_1)) );
            lin_2_sample = normrnd( lin_2, patch_linacc_std(sens(2))*sqrt(sum(lin_2.*lin_2)) );
            lin_3_sample = normrnd( lin_3, patch_linacc_std(sens(3))*sqrt(sum(lin_3.*lin_3)) );
        else
            lin_1_sample = normrnd( lin_1, patch_linacc_std(sens(1)) );
            lin_2_sample = normrnd( lin_2, patch_linacc_std(sens(2)) );
            lin_3_sample = normrnd( lin_3, patch_linacc_std(sens(3)) );
        end
        
        lin_a = lin_1_sample'; lin_aT = patch_transforms{sens(1)};
        lin_b = lin_2_sample'; lin_bT = patch_transforms{sens(2)};
        lin_c = lin_3_sample'; lin_cT = patch_transforms{sens(3)};
        
        particle_states(i,4:6) = ComputeAngAccel( ang_vel_sample', lin_a, lin_b, lin_c, lin_aT, lin_bT, lin_cT )';
        particle_states(i,7:9) = ProjectLinAcc( lin_a, ang_vel_sample, particle_states(i,4:6), lin_aT );
    end
    
    if ( std_scale == 1 )
        scaled_angvel_std = patch_angvel_std;
        for i=1:ndat
            mag_av = sqrt( sum( patch_data{i}.ang_vel(1,:) .* patch_data{i}.ang_vel(1,:), 2 ) );
            %mag_av = sqrt( sum( particle_states(i,1:3) .* particle_states(i,1:3) ) );
            scaled_angvel_std(i,:) = patch_angvel_std(i,:) * mag_av;
        end
        particle_weights = CalculateWeights( particle_states(:,1:3), all_angvels, scaled_angvel_std );
    else
        particle_weights = CalculateWeights( particle_states(:,1:3), all_angvels, patch_angvel_std );
    end
    
    % Check if need to resample particles
    n_eff = 1 / ( particle_weights' * particle_weights );
    if n_eff < n_particles / 2
        particle_states = ResampleParticles(particle_states, particle_weights);
    end
    
    state_mat(1,:) = particle_weights' * particle_states;
    for i=1:size( state_mat, 2 )
        state_std(1,i) = particle_weights' * ( particle_states(:,i) - state_mat(1,i) ).^2;
    end
    
%     if debug == 1
%         figure(1); clf;
%         curr_t = t_vec(1);
%         subplot(1,3,1); hold on;
%         plot( ones(n_particles,1)*curr_t, particle_states(:,1), 'bx' )
%         plot( [curr_t, curr_t], [ang_vel_filt(1,1) - ang_vel_std(1,1), ang_vel_filt(1,1) + ang_vel_std(1,1)], 'k-' );
%         plot( curr_t, ang_vel_filt(1,1), 'r+' )
%         xlim([t_vec(1), t_vec(end)])
%         
%         subplot(1,3,2); hold on;
%         plot( ones(n_particles,1)*curr_t, particle_states(:,2), 'bx' )
%         plot( [curr_t, curr_t], [ang_vel_filt(1,2) - ang_vel_std(1,2), ang_vel_filt(1,2) + ang_vel_std(1,2)], 'k-' );
%         plot( curr_t, ang_vel_filt(1,2), 'r+' )
%         xlim([t_vec(1), t_vec(end)])
%         
%         subplot(1,3,3); hold on;
%         plot( ones(n_particles,1)*curr_t, particle_states(:,3), 'bx' )
%         plot( [curr_t, curr_t], [ang_vel_filt(1,3) - ang_vel_std(1,3), ang_vel_filt(1,3) + ang_vel_std(1,3)], 'k-' );
%         plot( curr_t, ang_vel_filt(1,3), 'r+' )
%         xlim([t_vec(1), t_vec(end)])
%         keyboard;
%     end
    
    %% Particle filter to determine states
    for i=2:length( t_vec )
        
        % Estimate angular acceleration
        dt = t_vec(i) - t_vec(i-1);
        for j=1:n_particles
            % Propagate
            ang_vel_sample = particle_states(j,1:3);
            particle_states(j,1:3) = particle_states(j,4:6) * dt + ang_vel_sample;
            
            % Compute angular acceleration state estimate - Should this
            % happen after resampling?
            %sens = randperm(ndat);
            sens = circshift((1:ndat)',randi([0,ndat-1]));
            ang_vel_sample = particle_states(j,1:3);
        
            lin_1 = patch_data{sens(1)}.lin_acc(i,:);
            lin_2 = patch_data{sens(2)}.lin_acc(i,:);
            lin_3 = patch_data{sens(3)}.lin_acc(i,:);
            
            if ( std_scale == 1 )
                lin_1_sample = normrnd( lin_1, patch_linacc_std(sens(1))*sqrt(sum(lin_1.*lin_1)) );
                lin_2_sample = normrnd( lin_2, patch_linacc_std(sens(2))*sqrt(sum(lin_2.*lin_2)) );
                lin_3_sample = normrnd( lin_3, patch_linacc_std(sens(3))*sqrt(sum(lin_3.*lin_3)) );
            else
                lin_1_sample = normrnd( lin_1, patch_linacc_std(sens(1)) );
                lin_2_sample = normrnd( lin_2, patch_linacc_std(sens(2)) );
                lin_3_sample = normrnd( lin_3, patch_linacc_std(sens(3)) );
            end
        
            lin_a = lin_1_sample'; lin_aT = patch_transforms{sens(1)};
            lin_b = lin_2_sample'; lin_bT = patch_transforms{sens(2)};
            lin_c = lin_3_sample'; lin_cT = patch_transforms{sens(3)};
            
            particle_states(j,4:6) = ComputeAngAccel( ang_vel_sample', lin_a, lin_b, lin_c, lin_aT, lin_bT, lin_cT )';
            particle_states(j,7:9) = ProjectLinAcc( lin_a, ang_vel_sample, particle_states(j,4:6), lin_aT );
        end
        
        % Observations
        all_angvels = zeros(ndat,3);
        for j=1:ndat
            all_angvels(j,:) = GetCurrAngVel( patch_data{j}, patch_transforms{j}, i );
        end
        
        if ( std_scale == 1 )
            scaled_angvel_std = patch_angvel_std;
%             for j=1:ndat
%                mag_av = sqrt( sum( patch_data{j}.ang_vel(i,:) .* patch_data{j}.ang_vel(i,:), 2 ) );
%                %scaled_angvel_std(j,:) = patch_angvel_std(j,:) .* abs( patch_data{j}.ang_vel(i,:) )
%                %mag_av = sqrt( sum( state_mat(i-1,1:3) .* state_mat(i-1,1:3) ) );
%                scaled_angvel_std(j,:) = patch_angvel_std(j,:) * mag_av;
%             end
            mag_av = sqrt( sum( state_mat(i-1,1:3) .* state_mat(i-1,1:3) ) );
            
            if (mag_av > 0.5)
                scaled_angvel_std = patch_angvel_std * mag_av;
            end
            
            particle_weights = CalculateWeights( particle_states(:,1:3), all_angvels, scaled_angvel_std );
        else
            particle_weights = CalculateWeights( particle_states(:,1:3), all_angvels, patch_angvel_std );
        end

        % Check if need to resample Particles
        n_eff = 1 / ( particle_weights' * particle_weights );
        if n_eff < n_particles / 2
           particle_states = ResampleParticles(particle_states, particle_weights);
           if ( std_scale == 1 )
               particle_weights = CalculateWeights( particle_states(:,1:3), all_angvels, scaled_angvel_std );
           else
               particle_weights = CalculateWeights( particle_states(:,1:3), all_angvels, patch_angvel_std );
           end
        end
        
        % Estimate
        state_mat(i,:) = particle_weights' * particle_states;
        for j=1:size( state_mat, 2 )
            state_std(i,j) = particle_weights' * ( particle_states(:,j) - state_mat(i,j) ).^2;
        end
         
%         if debug == 1
%             figure(1); clf;
%             curr_t = t_vec(i);
%             subplot(1,3,1); hold on;
%             plot( ones(n_particles,1)*curr_t, particle_states(:,1), 'bx' )
%             plot( [curr_t, curr_t], [state_mat(i,1) - state_std(i,1), state_mat(i,1) + state_std(i,1)], 'k-' );
%             plot( t_vec(1:i), state_mat(1:i,1), 'r+' )
%             xlim([t_vec(1), t_vec(end)])
% 
%             subplot(1,3,2); hold on;
%             plot( ones(n_particles,1)*curr_t, particle_states(:,2), 'bx' )
%             plot( [curr_t, curr_t], [state_mat(i,2) - state_std(i,2), state_mat(i,2) + state_std(i,2)], 'k-' );
%             plot( t_vec(1:i), state_mat(1:i,2), 'r+' )
%             xlim([t_vec(1), t_vec(end)])
% 
%             subplot(1,3,3); hold on;
%             plot( ones(n_particles,1)*curr_t, particle_states(:,3), 'bx' )
%             plot( [curr_t, curr_t], [state_mat(i,3) - state_std(i,3), state_mat(i,3) + state_std(i,3)], 'k-' );
%             plot( t_vec(1:i), state_mat(1:i,3), 'r+' )
%             xlim([t_vec(1), t_vec(end)])
%             keyboard;
%        end
    end
    
    %% Based on velocity, recalculate other state variables
    % Angular acceleration is derivative of the angular velocity
    
    % Linear acceleration is a weighted average of all the measured linear
    % accelerations.
    
    filtered_data.t = t_vec;
    filtered_data.ang_vel = state_mat(:,1:3);
    filtered_data.ang_acc = state_mat(:,4:6);
    filtered_data.lin_acc = state_mat(:,7:9);
    filtered_data.ang_vel_std = state_std(:,1:3);
    filtered_data.ang_acc_std = state_std(:,4:6);
    filtered_data.lin_acc_std = state_std(:,7:9);
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