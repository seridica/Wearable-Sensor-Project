%% Header
% AA228 Final Project
% Calvin Kuo
% 5/11/2016
% ImprovedKinematicsParticleFilterNew.m

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

function filtered_data = ImprovedKinematicsParticleFilterNew( patch_data, patch_transforms, in_patch_stds, std_scale )
    
    % Scaling error
    if nargin == 3
        std_scale = 0;
    end

    % Flag for debugging
    debug = 0;

    assert( length( patch_data ) == length( patch_transforms ) );
    assert( length( patch_data ) == length( in_patch_stds ) );
    
    ndat = length( patch_data );
    
    % Pull useful values out
    patch_stds = ones(ndat,6);
    for i=1:ndat
        patch_stds(i,4:6) = in_patch_stds{i}(1,:);
        patch_stds(i,1:3) = in_patch_stds{i}(2,:);
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
    all_obs = zeros(ndat,6);
    for i=1:ndat
        all_obs(i,:) = GetCurrMeas( patch_data{i}, patch_transforms{i}, 1 );
    end
    
    % Initialize Particles
    n_particles = 10000;
    
    % Initialize angular velocity and linear acceleration estimates based
    % on observed values.
    % Angular acceleration give a uniform distribution over all three axes between -5000, 5000 rad/s^2
    % Line acceleration give a uniform distribution over all three axes
    % between -30, 30
    mean_angvels = mean( all_obs(:,1:3), 1 );
    std_angvels = std( all_obs(:,1:3), 1 );
    particle_states = zeros( n_particles, nstates );
    %for i=1:n_particles
    %    particle_states(i,1:3) = mean_angvels - std_angvels + 2*std_angvels.*rand(1,3);
%         particle_states(i,4:6) = -500+rand(3,1)*1000;
%         particle_states(i,7:9) = -5+rand(3,1)*10;
    %    particle_states(i,4:6) = -100+rand(3,1)*200;
    %    particle_states(i,7:9) = -1+rand(3,1)*2;
    %end
    particle_states(:,1:3) = repmat( mean_angvels, n_particles, 1 ) + repmat( std_angvels, n_particles, 1 ) .* ( -1 + 2*rand(n_particles,3) );
    particle_states(:,4:6) = -100+rand(n_particles, 3)*200;
    particle_states(:,7:9) = -1+rand(n_particles, 3)*2;
    
    if ( std_scale == 1 )
        %scaled_stds = patch_stds;
        scaled_stds = repmat( [2.0, 2.0, 2.0, 0.5, 0.5, 0.5], 4, 1 );
        particle_weights = CalculateWeightsNew( particle_states, all_obs, scaled_stds, patch_transforms );
    else
        particle_weights = CalculateWeightsNew( particle_states, all_obs, patch_stds, patch_transforms );
    end
    
    % Check if need to resample particles
    n_eff = 1 / ( particle_weights' * particle_weights );
    %keyboard
    if n_eff < n_particles / 2
        [particle_states, particle_weights]  = ResampleParticles(particle_states, particle_weights);
    end
    
    % EST
%     temp_est = particle_weights' * particle_states;
%     temp_std = temp_est;
%     for j=1:size( state_mat, 2 )
%         temp_std(j) = particle_weights' * ( particle_states(:,j) - temp_est(j) ).^2;
%     end

%         for j=1:3
%            %particle_states(:,j) = normrnd( temp_est(j), max(0.5, sqrt(temp_std(j))), [n_particles, 1] );
%            particle_states(:,j) = normrnd( temp_est(j), 0.1, [n_particles, 1] );
%         end
%         for j=4:6
%             %particle_states(:,j) = normrnd( temp_est(j), max(50, sqrt(temp_std(j))), [n_particles, 1] );
%             particle_states(:,j) = normrnd( temp_est(j), 50 , [n_particles, 1] );
%         end
%         for j=7:9
%             %particle_states(:,j) = normrnd( temp_est(j), max(0.5, sqrt(temp_std(j))), [n_particles, 1] );
%             particle_states(:,j) = normrnd( temp_est(j), 0.5, [n_particles, 1] );
%         end
% 
%     for j=1:3
%        particle_states(:,j) = particle_states(:,j) + normrnd( 0, 0.1, [n_particles, 1] );
%     end
%     for j=4:6
%         particle_states(:,j) = particle_states(:,j) + normrnd( 0, 10, [n_particles, 1] );
%     end
%     for j=7:9
%         particle_states(:,j) = particle_states(:,j) + normrnd( 0, 0.1, [n_particles, 1] );
%     end


%     if (std_scale == 1)
%         particle_weights = CalculateWeightsNew( particle_states, all_obs, scaled_stds, patch_transforms );
%     else
%         particle_weights = CalculateWeightsNew( particle_states, all_obs, patch_stds, patch_transforms );
%     end


    n_eff = 1 / ( particle_weights' * particle_weights );
    %keyboard;
    if n_eff < n_particles / 2
       [particle_states, particle_weights] = ResampleParticles(particle_states, particle_weights);
    end
    
    
    state_mat(1,:) = particle_weights' * particle_states;
    disp( length( unique( particle_weights ) ) );
    disp( state_mat(i,:) );
    for i=1:size( state_mat, 2 )
        state_std(1,i) = particle_weights' * ( particle_states(:,i) - state_mat(1,i) ).^2;
    end
    %keyboard;
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
        
        mag_aa = sqrt( sum( state_mat(i-1,4:6) .* state_mat(i-1,4:6) ) );
        mag_la = sqrt( sum( state_mat(i-1,7:9) .* state_mat(i-1,7:9) ) );
        
        % Estimate angular acceleration
        dt = t_vec(i) - t_vec(i-1);
%         for j=1:n_particles
%             % Propagate
%             ang_vel_sample = particle_states(j,1:3);
%             particle_states(j,1:3) = particle_states(j,4:6) * dt + ang_vel_sample;
%             % Resample
%             
%             particle_states(j,4:6) = state_mat(i-1,4:6)'-[200;200;200]+rand(3,1)*400;
%             particle_states(j,7:9) = state_mat(i-1,7:9)'-[2;2;2]+rand(3,1)*4;
% %             particle_states(j,4:6) = -500+rand(3,1)*1000;
% %             particle_states(j,7:9) = -5+rand(3,1)*10;
%         end

%         particle_states(:,1:3) = particle_states(:,1:3) + dt*particle_states(:,4:6) - 0.1 + 0.2*rand(n_particles,3) ;
%         particle_states(:,4:6) = repmat( state_mat(i-1,4:6), n_particles, 1 ) - 500 + 1000*rand(n_particles,3);
%         particle_states(:,7:9) = repmat( state_mat(i-1,7:9), n_particles, 1 ) - 5 + 10*rand(n_particles,3);

        %particle_states(:,1:3) = particle_states(:,1:3) + dt*particle_states(:,4:6) + normrnd(0, 0.2, [n_particles,3]) ;
        %particle_states(:,4:6) = repmat( state_mat(i-1,4:6), n_particles, 1 ) + normrnd(0, 200, [n_particles,3]);
        %particle_states(:,7:9) = repmat( state_mat(i-1,7:9), n_particles, 1 ) + normrnd(0, 2, [n_particles,3]);

        particle_states(:,1:3) = particle_states(:,1:3) + dt*particle_states(:,4:6) + normrnd(0, 0.1, [n_particles,3]) ;
        particle_states(:,4:6) = repmat( state_mat(i-1,4:6), n_particles, 1 ) + normrnd(0, 20, [n_particles,3]);
        particle_states(:,7:9) = repmat( state_mat(i-1,7:9), n_particles, 1 ) + normrnd(0, 1, [n_particles,3]);
        
        % Observations
        all_obs = zeros(ndat,6);
        for j=1:ndat
            all_obs(j,:) = GetCurrMeas( patch_data{j}, patch_transforms{j}, i );
        end
        
        if ( std_scale == 1 )
%             for j=1:ndat
%                mag_av = sqrt( sum( patch_data{j}.ang_vel(i,:) .* patch_data{j}.ang_vel(i,:), 2 ) );
%                %scaled_angvel_std(j,:) = patch_angvel_std(j,:) .* abs( patch_data{j}.ang_vel(i,:) )
%                %mag_av = sqrt( sum( state_mat(i-1,1:3) .* state_mat(i-1,1:3) ) );
%                scaled_angvel_std(j,:) = patch_angvel_std(j,:) * mag_av;
%             end
            mag_av = sqrt( sum( state_mat(i-1,1:3) .* state_mat(i-1,1:3) ) );
            
            scaled_stds(:,1:3) = max( patch_stds(:,1:3) * mag_av, 2.0 );
            scaled_stds(:,4:6) = max( patch_stds(:,4:6) * mag_la, 0.5 );
            
            particle_weights = CalculateWeightsNew( particle_states, all_obs, scaled_stds, patch_transforms );
        else
            particle_weights = CalculateWeightsNew( particle_states, all_obs, patch_stds, patch_transforms );
        end

        % Check if need to resample Particles
%         n_eff = 1 / ( particle_weights' * particle_weights );
%         %keyboard;
%         if n_eff < n_particles / 2
%            [particle_states, particle_weights] = ResampleParticles(particle_states, particle_weights);
%         end
        
        % Estimate
%         temp_est = particle_weights' * particle_states;
%         temp_std = temp_est;
%         for j=1:size( state_mat, 2 )
%             temp_std(j) = particle_weights' * ( particle_states(:,j) - temp_est(j) ).^2;
%         end
% 
%         for j=1:3
%             particle_states(:,j) = normrnd( temp_est(j), max(0.5, sqrt(temp_std(j))), [n_particles, 1] );
%             %disp(max(0.5, sqrt(temp_std(j))));
%             %particle_states(:,j) = normrnd( temp_est(j), 0.1, [n_particles, 1] );
%         end
%         for j=4:6
%             particle_states(:,j) = normrnd( temp_est(j), max(50, sqrt(temp_std(j))), [n_particles, 1] );
%             %disp(max(50, sqrt(temp_std(j))));
%             %particle_states(:,j) = normrnd( temp_est(j), 50 , [n_particles, 1] );
%         end
%         for j=7:9
%             particle_states(:,j) = normrnd( temp_est(j), max(0.5, sqrt(temp_std(j))), [n_particles, 1] );
%             %disp(max(0.5, sqrt(temp_std(j))));
%             %particle_states(:,j) = normrnd( temp_est(j), 0.5, [n_particles, 1] );
%         end
% 
%         for j=1:3
%            particle_states(:,j) = particle_states(:,j) + normrnd( 0, 0.1, [n_particles, 1] );
%         end
%         for j=4:6
%             particle_states(:,j) = particle_states(:,j) + normrnd( 0, 10, [n_particles, 1] );
%         end
%         for j=7:9
%             particle_states(:,j) = particle_states(:,j) + normrnd( 0, 0.1, [n_particles, 1] );
%         end


%         if (std_scale == 1)
%             particle_weights = CalculateWeightsNew( particle_states, all_obs, scaled_stds, patch_transforms );
%         else
%             particle_weights = CalculateWeightsNew( particle_states, all_obs, patch_stds, patch_transforms );
%         end
        
        
        n_eff = 1 / ( particle_weights' * particle_weights );
        %keyboard;
        if n_eff < n_particles / 2
           [particle_states, particle_weights] = ResampleParticles(particle_states, particle_weights);
        end
        
        
        state_mat(i,:) = particle_weights' * particle_states;
        
        disp( length( unique( particle_weights ) ) );
        disp( state_mat(i,:) );
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
%        keyboard;
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

%% Helper function that obtains the current measurement vector
% measurement vector at time t
function curr_meas = GetCurrMeas( indata, inT, curr_t )
    this_angvel = indata.ang_vel(curr_t,:);
    this_linacc = indata.lin_acc(curr_t,:);
    R = inT(1:3,1:3);
    curr_angvel = this_angvel*R;
    curr_linacc = this_linacc*R;
    
    curr_meas = [curr_angvel, curr_linacc];
end

%% Helper function that resamples particles based on distribution
% Stratified Resampling
% https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python/blob/master/12-Particle-Filters.ipynb
function [new_particles, new_weights] = ResampleParticles( curr_particles, curr_weights )
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
    new_weights = curr_weights(indexes,:);
    new_weights = new_weights / ( sum( new_weights ) );
end

% Systematic Resampling
% https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python/blob/master/12-Particle-Filters.ipynb