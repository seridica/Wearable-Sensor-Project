%% Header
% File name: StitchUBC.m
% Author: Calvin Kuo
% Date: 10/03/2018

% Function stitches together data from a continuous trial according to the
% UBC firmware.

function processed_data = StitchUBC( impacts )
    
    % Set up data structures
    t = [];
    lin_acc = [];
    ang_vel = [];

    % Loop through remaining impacts
    for i=1:( length( impacts ) - 1)
        
        % Get timings
        t1 = impacts(i).Info.t1;
        t2 = impacts(i).Info.t2;
        t3 = impacts(i).Info.t3;
        t4 = impacts(i).Info.t4;
        
        % t will start however many milliseconds after the impact (t2-t1)
        if ( i == 1 )
            t_begin_rst = t3 - t1;
            t_begin_hsi = double( t4 - t2 );
            t_init_rst = t3;
            t_init_hsi = double( t4 );
            lin_acc = [lin_acc; impacts(i).lin_acc(43:84,:)];
            
            temp_av = impacts(i).ang_vel(43:84,:);
            temp_av(find(temp_av(:,1) > 40),:) = [];
            ang_vel = [ang_vel; temp_av];
        elseif (i == length(impacts)-1)
            lin_acc = [lin_acc; impacts(i).lin_acc(1:21,:)];
            
            temp_av = impacts(i).ang_vel(1:21,:);
            temp_av(find(temp_av(:,1) > 40),:) = [];
            ang_vel = [ang_vel; temp_av];
            t_end_rst = t1;
            t_end_hsi = double( t2 );
        else
            lin_acc = [lin_acc; impacts(i).lin_acc(1:84,:)];
            
            temp_av = impacts(i).ang_vel(1:84,:);
            temp_av(find(temp_av(:,1) > 40),:) = [];
            ang_vel = [ang_vel; temp_av];
        end
    end
    
    % Make time vector
    nSamples = size( lin_acc, 1 );
    dt = (t_end_rst - t_init_rst) / nSamples / 1000;
    t_begin_rst = t_begin_rst / 1000;
    t = t_begin_rst:dt:(t_begin_rst+(nSamples-1)*dt);
    
    avSamples = size( ang_vel, 1 );
    dtav = (t_end_rst - t_init_rst) / avSamples / 1000;
    tav = t_begin_rst:dtav:(t_begin_rst+(avSamples-1)*dtav);
    ang_vel = interp1( tav, ang_vel, t, 'spline' );
    
    % Filter at 300Hz for lin acc
    fs = (1/dt);
    [b,a] = butter(2, 150 / ( fs/2 ));
    filt_linacc = filtfilt( b,a, lin_acc );
    
    % Filter at 110Hz for ang vel
    fs = (1/dt);
    [b,a] = butter(2, 150 / ( fs/2 ));
    filt_angvel = filtfilt( b,a, ang_vel );
    
    % Output
    processed_data.t = t;
    processed_data.lin_acc = filt_linacc;
    processed_data.ang_vel = filt_angvel;
end