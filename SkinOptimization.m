%%%
% File name: SkinOptimization.m
% Author: Calvin Kuo
% Date: 10/30/2018
%
% Optimizes for various skin models

function [optMotion, optParams] = SkinOptimization( t_in, motion_in, motion_comp, init_params, model_type )
    switch model_type
        case 1
            lb = [0.1, 0.01, 0.005];
            ub = [5000.0, 50.0, 1.0];
            opts = optimoptions('fmincon','Display','iter');
            optParams = fmincon(@(x) LinearSpringDamperOpt(t_in, motion_in, motion_comp, x), init_params, [], [], [], [], lb, ub, [], opts);
        case 2
            lb = []
    end
    
    optMotion = zeros( size( motion_in ) );
    for j=1:size( motion_in, 2 )
        optMotion(:,j) = SkinModel( t_in, motion_in(:,j), optParams, model_type );

        figure(j); clf; hold on;
        plot( t_in, motion_in(:,j) );
        plot( t_in, motion_comp(:,j) );
        plot( t_in, optMotion(:,j) );
    end
end

function cost = LinearSpringDamperOpt( t_in, motion_in, motion_comp, param_in )
    
    nImpacts = size( motion_in, 2 );
    cost = 0;
    for i=1:nImpacts
        outMotion = SkinModel( t_in, motion_in(:,i), param_in, 1 );
        cost = cost + norm( outMotion - motion_comp(:,i) );
    end
    
end