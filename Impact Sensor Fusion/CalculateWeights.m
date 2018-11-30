%% Header
% AA228 Final Project
% Calvin Kuo
% 11/27/2016
% CalculateWeights.m

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
function particle_weights = CalculateWeights( curr_particles, curr_observations, observation_stds )

    % Check on size
    assert( size( curr_particles, 2 ) == size( curr_observations, 2 ) );
    
    % Determine weights of the current particles based on observations
    nparticles = size( curr_particles, 1 );
    nobservations = size( curr_observations, 1 );
    nstates = size( curr_particles, 2 );
    
    % Make probability distributions for the observations
    clear('obs_pdfs');
    for i=1:nobservations
        for j=1:nstates
            %observation_stds
            %curr_observations
            obs_pdfs(i,j) = makedist( 'Normal', curr_observations(i,j), observation_stds(i,j) );
        end
    end
    
    % Compute weights based on probability distributions
    particle_weights = ones( nparticles, 1 );
    for i=1:nobservations
        for j=1:nstates
            this_pdf = obs_pdfs(i,j);
            particle_weights = particle_weights .* pdf( this_pdf, curr_particles(:,j) );
        end
    end
    particle_weights = particle_weights / ( sum( particle_weights ) );
    
end