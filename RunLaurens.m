%% Helper function for Laurens Model
function [estimate_velocity, stored_velocity] = RunLaurens( t, input_angvel, input_angacc, input_linacc )
    
    %% PARAMETERS
    params.tc = 7.3; %4;
    params.tn = 75;
    params.kv = 0.52; %0.247;
    params.tvs = 7.7; %14;

    % Switches
    params.vswitch = 0;
    params.gswitch = 0;

    params.kf = 0.0; % Set to zero to turn off gravity
    params.ts = 1/0.65;

    params.go = 0; % Set to zero to turn off vision
    params.ko = 0; % Set to zero to turn off vision

    params.gd = 1;
    
    %% INPUTS
    u.alpha = [t'; input_angacc'];
    u.omega = [t'; input_angvel'];
    u.acc = [t'; input_linacc' * 9.81];

    % Gravity
    u.grav = [0; -9.81; 0];

    % Angular position (for now)
    u.Tcan = [1 0 0;...
              0 1 0;...
              0 0 1];
          
    u.Cnoise = [0;0;0];
    u.dVSnoise = [0;0;0];
    u.VInoise = [0;0;0];
    
    %% Initializing initial conditions
    C = input_angvel(1,:)';              % Canal output signal   
    D = input_angvel(1,:)';              % Cana afferent signal
    INT = [0;0;0];            % Integral output (Endolymph)
    VS = [0;0;0];             % Velocity storage leaky integrator
    VSf = [0;0;0];            % Velocity storage full
    GE = u.grav;             % Gravity estimate

    init_x = [C; D; INT; VS; VSf; GE];

    [tinteg, full_states] = ode45(@(t,y) LaurensVestibularModelDeriv(t,y,u,params), t, init_x);
    % Get velocity estimate
    nSteps = length( tinteg );
    estimate_velocity = zeros( nSteps, 4 );
    estimate_velocity(:,1) = tinteg;
    for i=1:length(tinteg)
        u.D = full_states(i,4:6)';
        u.dirVI = [0;0;0];
        u.VSf = full_states(i,13:15)';
        estimate_velocity(i,2:4) = VelocityEstimate( tinteg(i), u, params )';
    end
    stored_velocity = [tinteg, full_states(:,10:12)];
end