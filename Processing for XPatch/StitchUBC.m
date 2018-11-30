%% Header
% File name: StitchUBC.m
% Author: Calvin Kuo
% Date: 10/03/2018

% Function stitches together data from a continuous trial according to the
% UBC firmware.

function processed_data = StitchUBC( impacts, startImpact, offset, scale_angvel )
    
    % Set up data structures
    t = [];
    lin_acc = [];
    ang_vel = [];
    
    % Get timings
    t0p = impacts(1).Info.t2;
    t0 = t0p / 1000;
    t1 = impacts(startImpact).Info.t1;
    t2 = impacts(startImpact).Info.t2;
    t3 = impacts(startImpact).Info.t3;
    t4 = double( impacts(startImpact).Info.t4 );
    if (offset == 1)
        t_begin_rst = ( t2 - t0p ) / 1000;
        t_init_rst = t2;
        lin_acc = [lin_acc; impacts(startImpact).lin_acc(43:84,:)];
            
        temp_av = impacts(startImpact).ang_vel(43:84,:);
        temp_av(find(temp_av(:,1) > 40),:) = temp_av(find(temp_av(:,1) > 40)-1,:);
        ang_vel = [ang_vel; temp_av * scale_angvel];
        
        initInd = startImpact + 1;
    elseif (offset == 2)
        t_begin_rst = ( t3 - t0p ) / 1000;
        t_init_rst = t3;
        lin_acc = [lin_acc; impacts(startImpact).lin_acc(64:84,:)];
            
        temp_av = impacts(startImpact).ang_vel(64:84,:);
        temp_av(find(temp_av(:,1) > 40),:) = temp_av(find(temp_av(:,1) > 40)-1,:);
        ang_vel = [ang_vel; temp_av * scale_angvel];
        
        initInd = startImpact + 1;
        
        tpiece = linspace( t3, t4, 21 ) / 1000;
        t = [ t, tpiece - t0 ];
    elseif (offset == 3)
        t_begin_rst = ( t4 - t0p ) / 1000;
        t_init_rst = t4;
        
        initInd = startImpact + 1;
    else
        t_begin_rst = ( impacts(startImpact+1).Info.t1 - t0p ) / 1000;
        t_init_rst = impacts(startImpact+1).Info.t1;
        lin_acc = [lin_acc; impacts(startImpact+1).lin_acc(22:84,:)];
        
        temp_av = impacts(startImpact+1).ang_vel(22:84,:);
        temp_av(find(temp_av(:,1) > 40),:) = temp_av(find(temp_av(:,1) > 40)-1,:);
        ang_vel = [ang_vel; temp_av * scale_angvel];
        
        initInd = startImpact + 2;
    end
    t_prev = t4;

    % Loop through remaining impacts
    for i=initInd:( length( impacts ) - 13)
        
        % Get timings
        t1 = impacts(i).Info.t1;
        t2 = impacts(i).Info.t2;
        t3 = impacts(i).Info.t3;
        t4 = double( impacts(i).Info.t4 );
        
        if (i == length(impacts)-13)
            lin_acc = [lin_acc; impacts(i).lin_acc(1:84,:)];
            
            temp_av = impacts(i).ang_vel(1:84,:);
            temp_av(find(temp_av(:,1) > 40),:) = temp_av(find(temp_av(:,1) > 40)-1,:);
            ang_vel = [ang_vel; temp_av * scale_angvel];
            t_end_rst = double( t4 )
        else
            lin_acc = [lin_acc; impacts(i).lin_acc(1:84,:)];
            
            temp_av = impacts(i).ang_vel(1:84,:);
            temp_av(find(temp_av(:,1) > 40),:) = temp_av(find(temp_av(:,1) > 40)-1,:);
            ang_vel = [ang_vel; temp_av * scale_angvel];
        end
        
        tpiece = linspace( t_prev, t1, 22 ) / 1000;
        t = [ t, tpiece(2:end) - t0 ];
        
        tpiece = linspace( t1, t2, 22 ) / 1000;
        t = [ t, tpiece(2:end) - t0 ];
        
        tpiece = linspace( t2, t3, 22 ) / 1000;
        t = [ t, tpiece(2:end) - t0 ];
        
        tpiece = linspace( t3, t4, 22 ) / 1000;
        t = [ t, tpiece(2:end) - t0 ];
        
        t_prev = t4;
    end
    
    % Make time vector
    nSamples = size( lin_acc, 1 );
    dt = (t_end_rst - t_init_rst) / nSamples / 1000
    t_begin_rst = t_begin_rst / 1000;
    %t = t_begin_rst:dt:(t_begin_rst+(nSamples-1)*dt);
    
    %avSamples = size( ang_vel, 1 );
    %dtav = (t_end_rst - t_init_rst) / avSamples / 1000;
    %tav = t_begin_rst:dtav:(t_begin_rst+(avSamples-1)*dtav);
    %ang_vel = interp1( tav, ang_vel, t, 'spline' );
    
    % Filter at 300Hz for lin acc
    fs = (1/dt)
    [b,a] = butter(2, 150 / ( fs/2 ));
    filt_linacc = filtfilt( b,a, lin_acc );
    
    % Filter at 110Hz for ang vel
    fs = (1/dt);
    [b,a] = butter(2, 110 / ( fs/2 ));
    filt_angvel = filtfilt( b,a, ang_vel );
    
    % Output
    processed_data.t = t;
    processed_data.lin_acc = filt_linacc;
    processed_data.ang_vel = filt_angvel;
end