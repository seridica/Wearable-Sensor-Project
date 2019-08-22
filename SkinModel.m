%%%
% File name: SkinModels.m
% Author: Calvin Kuo
% Date: 10/30/2018
%
% Various skin models

function outMotion = SkinModel( t_in, motion_in, params, model_type )
    switch model_type
        case 1
            outMotion = LinearSpringDamper( t_in, motion_in, params );
        case 2
            outMotion = PendulumCartSpringDamper( t_in, motion_in, params );
        case 3
            outMotion = LinearSpringDamperToe( t_in, motion_in, params );
    end
end

function outMotion = LinearSpringDamper( t_in, motion_in, params )
    k = params(1);
    b = params(2);
    m = params(3);
    
    u = [t_in, motion_in*9.81];
    
    opts = odeset('MaxStep', 0.0025);
    [tmotion, full_states] = ode45(@(t,y) LinearSpringDamperDeriv( t, y, k, b, m, u ), t_in, [0,0], opts);
    
    ddy = zeros( length( tmotion ), 1 );
    for i=1:length( tmotion )
        this_ddy = LinearSpringDamperDeriv( tmotion(i), full_states(i,:), k, b, m, u );
        ddy(i) = this_ddy(1);
    end
    %outMotion = ddy / 9.81;
    outMotion = ( motion_in*9.81 - ddy ) / 9.81;
end

function outMotion = PendulumCartSpringDamper( t_in, motion_in, params )
    k = params(1);
    b = params(2);
    m = params(3);
    l = params(4);
    
    u = [t_in, motion_in(:,1)*9.81, motion_in(:,2)];
    
    opts = odeset('MaxStep', 0.0025);
    [tmotion, full_states] = ode45(@(t,y) PendulumCartSpringDamperDeriv( t, y, k, b, m, l, u ), t_in, [0,0], opts);
    
    ddy = zeros( length( tmotion ), 1 );
    for i=1:length( tmotion )
        this_ddy = PendulumCartSpringDamperDeriv( tmotion(i), full_states(i,:), k, b, m, l, u );
        ddy(i) = this_ddy(1);
    end
    ddth = motion_in(:,2) - ddy;
    outMotion = [ motion_in(:,1) - ( l*ddth ) / 9.81, ddth, full_states(:,1) ];
    
    %outMotion = [ - ( l*ddth ) / 9.81, ddth ];
end

function outMotion = LinearSpringDamperToe( t_in, motion_in, params )
    k = params(1);
    b = params(2);
    m = params(3);
    toe = params(4);
    
    u = [t_in, motion_in*9.81];
    
    opts = odeset('MaxStep', 0.0025);
    [tmotion, full_states] = ode45(@(t,y) LinearSpringDamperToeDeriv( t, y, k, b, m, toe, u ), t_in, [0,0], opts);
    
    ddy = zeros( length( tmotion ), 1 );
    for i=1:length( tmotion )
        this_ddy = LinearSpringDamperToeDeriv( tmotion(i), full_states(i,:), k, b, m, toe, u );
        ddy(i) = this_ddy(1);
    end
    %outMotion = ddy / 9.81;
    outMotion = ( motion_in*9.81 - ddy ) / 9.81;
end

function dy = LinearSpringDamperDeriv( t, y, k, b, m, input_motion )
    u_in = interp1( input_motion(:,1), input_motion(:,2), t, 'spline' );

    dy = zeros(2,1);
    dy(2) = y(1);
    dy(1) = -k/m*y(2) - b/m*y(1) + u_in;
end

function dy = PendulumCartSpringDamperDeriv( t, y, k, b, m, l, input_motion )
    lin_in = interp1( input_motion(:,1), input_motion(:,2), t, 'spline' ); % Linear accel input
    ang_in = interp1( input_motion(:,1), input_motion(:,3), t, 'spline' ); % Angular accel input
    
    dy = zeros(2,1);
    dy(2) = y(1);
    dy(1) = -k/(m*l^2)*y(2) - b/(m*l^2)*y(1) + ang_in - lin_in / l;
end

function dy = LinearSpringDamperToeDeriv( t, y, k, b, m, d, input_motion )
    u_in = interp1( input_motion(:,1), input_motion(:,2), t, 'spline' );
    
    dy = zeros(2,1);
    dy(2) = y(1);
    if abs( y(2) ) < d
        dy(1) = -b/m*y(1) + u_in;
    elseif y(2) < -d
        dy(1) = -k/m*( y(2) + d )- b/m*y(1) + u_in;% * (b*abs(y(1))) + u_in; %
    else
        dy(1) = -k/m*( y(2) - d )- b/m*y(1) + u_in;% * (b*abs(y(1))) + u_in; %
    end
end