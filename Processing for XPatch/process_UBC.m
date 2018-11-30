function impact = process_UBC()

    % USER SELECTION (DEFAULT)
    % =========================================================================
    % Where is data file
    defaultDirectory = 'C:\Data\';
    [infiles, inpath] = uigetfile({'*.dat', 'Sensor data file (*.dat)'}, 'Select DAT file(s) with MG Sensor Data (OZ3)', defaultDirectory,'MultiSelect','on');

    % FILE IMPORT
    % =========================================================================    
    infiles = cellstr(infiles); % Convert infiles to cell array (to catch case of a single file)
    eventCount = 0; % Initialize event counter
    gravity = 9.80665;  % [m/s^2] Standard gravity
    
    % Cycle through each selected file
    for curFile = 1:length(infiles)

        % Open file
        pathAndFileName = [inpath infiles{curFile}];
        fid = fopen(pathAndFileName);
        if fid == -1, impact = []; break; end; % Break if fid is invalid
        nextLine = fgetl(fid); % Retrieve first line
        
        Info.File = infiles{curFile};

        % Reads impacts using readImpact function
        while ischar(nextLine)
            
            % DEVICE INFORMATION
            % =====================================================================
            % Update event counter
            eventCount = eventCount + 1;
            
            % Extra event information
            current_impact = readImpact_UBC(fid, nextLine)
            
            % Date information
            Info.ImpactNum = current_impact.impactnum
            
            % SAMPLING RATES
            % =====================================================================
            % Accelerometer sampling (to which gyroscope is upsampled)
            Info.SampleRate = 1000; % [Hz] Approximate
            Info.t1 = current_impact.t1;
            Info.t2 = current_impact.t2;
            Info.t3 = current_impact.t3;
            Info.t4 = current_impact.t4;
            Data.samples = 1:84; % [samples] Always the same number of samples: 84 (based on firmware)
            dt = 1/Info.SampleRate;
            Data.t = 0:dt:(length(Data.samples)*dt); % [s] Time vector, with 0 = time of trigger          
                
            % 1. LOAD RAW DATA
            % ===================================================================== 
            % Clear variables
            RawData = [];
            ShiftedData = [];
            TrimmedData = [];
            UpsampledData = [];
            OffsetData = [];
            ScaledData = [];
            TransformedData = [];
            FilteredData = [];
            ProjectedData = [];            
            
            % Raw data converted
            accelToG = 0.098 / 2; % 0.098 for 200g range, 0.098/2 for 100g range
            gyroToDegPerSec = 0.07; % V3
            Data.lin_acc = current_impact.lin_acceleration * accelToG;
            Data.ang_vel = current_impact.rot_velocity * gyroToDegPerSec * pi/180;
            
            impact(eventCount).Info = Info;
            impact(eventCount).Data = Data;
            impact(eventCount).t = Data.t';
            impact(eventCount).lin_acc = Data.lin_acc;         
            impact(eventCount).ang_vel = Data.ang_vel;

            % NEXT IMPACT
            % =====================================================================
            nextLine = fgetl(fid);  
        end
        
        % Close open files
        fclose('all');
    end
end

function output = readImpact_UBC(fid, line)

    % Initialize variables
    nextLine = line; 
    impactNum = nan; impactTime = nan;
    tick_impact = nan; tick_impactEnd = nan;
    gyroNum = nan; clack = nan; ims = nan; als = nan;
    linAcc_raw=[]; angVel_raw=[];

    % Parse the data file
    while ~strcmp(nextLine, '-----')
        category = textscan(nextLine, '%[^:]', 1);
        header = category{1}{1};
        switch header
            case 'Impact_number',          impactNum = cell2mat(textscan(nextLine,'Impact_number: %f'));
            case 'Event_time',             impactTime = cell2mat(textscan(nextLine,'Event_time: %f'));
            case 'Sys_tick_at_impact',     tick_impact = cell2mat(textscan(nextLine,'Sys_tick_at_impact: %f'));
            case 'Sys_tick_at_impact_end', tick_impactEnd = cell2mat(textscan(nextLine,'Sys_tick_at_impact_end: %f'));
            case 'Gyro_samples',           gyroNum = cell2mat(textscan(nextLine,'Gyro_samples: %f'));
            case 'Is_clack',               clack = cell2mat(textscan(nextLine,'Is_clack: %f'));
            case 'In_mouth_at_impact',     ims = cell2mat(textscan(nextLine,'In_mouth_at_impact: %f'));
            case 'ALS_at_impact',          als = cell2mat(textscan(nextLine,'ALS_at_impact: %f'));
            case 'Peak_linear_accel',      peak_linAcc = cell2mat(textscan(nextLine,'Peak_linear_accel: %f'));
            case 'Peak_rotational_accel',  peak_angAcc = cell2mat(textscan(nextLine,'Peak_rotational_accel: %f'));
            case 'CRC16',                  crc16 = cell2mat(textscan(nextLine,'CRC16: %f'));
            case 'Linear_data',            linAcc_raw = cell2mat(textscan(fgetl(fid),'[ %f, %f, %f ], '));
            case 'Rotational_data',        angVel_raw = cell2mat(textscan(fgetl(fid),'[ %f, %f, %f ], '));
        end
        nextLine=fgetl(fid);
    end
    
    % Calculate the magnitudes
    linAcc_mag_raw = sqrt(sum(linAcc_raw.^2,2));
    angVel_mag_raw = sqrt(sum(angVel_raw.^2,2));
    
    % Convert the time into useable format
    %[dateString, dateSerial] = x2time(impactTime);
    
    output.stringtime = '';
    output.serialtime = '';
    output.tick_impact = tick_impact;
    output.tick_impactEnd = tick_impactEnd; 
    output.lin_magnitude = linAcc_mag_raw;
    output.rot_magnitude = angVel_mag_raw;
    output.lin_acceleration = linAcc_raw;
    output.rot_velocity = angVel_raw;
    output.impactnum = impactNum;
    output.devicetime = impactTime;
    output.gyrosamples = gyroNum;
    
    output.t1 = impactTime;
    output.t2 = tick_impact;
    output.t3 = tick_impactEnd;
    output.t4 = bitshift( uint32(peak_angAcc), 16 ) + peak_linAcc;
    
    % HIJACKED output.isClack = crc16;
    output.gyrosamples_postimpact = crc16;
    
    % HIJACKED output.peak_linacc = peak_linAcc;
    output.impactindex_accel = peak_linAcc;
    
    % HIJACKED output.peak_angacc = peak_angAcc;
    output.impactindex_gyro = peak_angAcc;
    
    % HIJACKED output.hasALS = als;
    output.ir = als;   
    output.ims = ims;
end