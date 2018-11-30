%% Header
% Sensor Fusion Project
% Calvin Kuo
% 7/12/2017
% CalculateWeightsNew.m

%% File Description
% This file compares particles against observations to weight the particles

%% Function
% Inputs
% 1) Current particle states
% 2) Current observations of the states
% 3) Standard deviation for the observations
%
% Outputs
% 2) Particle weights

%% Helper function that 
function particle_weights = CalculateWeightsNew( curr_particles, curr_observations, observation_stds, patch_transforms )

    % Check on size
    assert( size( curr_observations, 2 ) == 6 );
    assert( size( observation_stds, 2 ) == 6 );
    assert( size( curr_observations, 1 ) == size( observation_stds, 1 ) );
    assert( size( curr_observations, 1 ) == length( patch_transforms ) );
    
    % Gravity
    g = 9.81;
    
    % Determine weights of the current particles based on observations
    nparticles = size( curr_particles, 1 );
    nobservations = size( curr_observations, 1 );
    nstates = size( curr_particles, 2 );
    
    % Make probability distributions for the observations
    clear('obs_pdfs');
    for i=1:nobservations
        %observation_stds
        %curr_observations
        % Observation Error
        obs_pdfs(i,1) = makedist( 'Normal', 0, observation_stds(i,1) );
        obs_pdfs(i,2) = makedist( 'Normal', 0, observation_stds(i,4) );
    end
    
    % Compute weights based on probability distributions
    particle_weights = ones( nparticles, 1 ); % * 1e300;
    for i=1:nparticles
        pang_vel = curr_particles(i,1:3)';
        pang_acc = curr_particles(i,4:6)';
        plin_acc = curr_particles(i,7:9)';
        for j=1:nobservations
            r_vec = patch_transforms{j}(1:3,4) * 2.54 / 100;
            
            % Project the particle linear acceleration to the observation
            % linear acceleration
            proj_acc = ( plin_acc*g + ( cross( pang_vel, cross( pang_vel, r_vec ) ) + cross( pang_acc, r_vec ) ) ) / g;
            %proj_acc = ( plin_acc*g + cross( pang_acc, r_vec ) ) / g;
            
%             this_pdf = obs_pdfs(j,1);
%             particle_weights(i) = particle_weights(i) * pdf( this_pdf, pang_vel(1) );
%             %particle_weights(i)
%             this_pdf = obs_pdfs(j,2);
%             particle_weights(i) = particle_weights(i) * pdf( this_pdf, pang_vel(2) );
%             %particle_weights(i)
%             this_pdf = obs_pdfs(j,3);
%             particle_weights(i) = particle_weights(i) * pdf( this_pdf, pang_vel(3) );
%             %particle_weights(i)
%             this_pdf = obs_pdfs(j,4);
%             particle_weights(i) = particle_weights(i) * pdf( this_pdf, proj_acc(1) );
%             %particle_weights(i)
%             this_pdf = obs_pdfs(j,5);
%             particle_weights(i) = particle_weights(i) * pdf( this_pdf, proj_acc(2) );
%             %particle_weights(i)
%             this_pdf = obs_pdfs(j,6);
%             particle_weights(i) = particle_weights(i) * pdf( this_pdf, proj_acc(3) );
%             %particle_weights(i)
% %             if ( particle_weights(i) < 1e-10 )
% %                 keyboard;
% %             end
            linacc_norm = norm( proj_acc - curr_observations(j,4:6)', 2 );
            angvel_norm = norm( pang_vel - curr_observations(j,1:3)', 2 );
            
            particle_weights(i) = particle_weights(i) * pdf( obs_pdfs(j,1), angvel_norm );
            particle_weights(i) = particle_weights(i) * pdf( obs_pdfs(j,2), linacc_norm );
% 
%             linacc_measurementErrorNorm = norm( proj_acc - curr_observations(j,4:6), 2 );
%             likelihood = 1/sqrt((2*pi).^3 * det(measurementNoise)) * exp(-0.5 * linacc_measurementErrorNorm);
%             particle_weights(i) = particle_weights(i) * likelihood;
        end
    end
    particle_weights = particle_weights / ( sum( particle_weights ) );
    
end