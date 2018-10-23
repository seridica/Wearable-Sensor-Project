function [dateString, dateSerial] = x2time(eventTime)
% [dateString, dateSerial] = x2time(eventTime)
% Computes serial and string formats of event time based on event time
% recorded on device (time since December 31, 2009 at 4 PM). The algorithm
% adjusts for Daylight Savings Time (DST) by adding an hour to the result
% if it occurred on a date considered DST.

    % Reference date
    x2start = datenum(2009,12,31,16,0,0);
    dateSerial = addtodate(x2start,eventTime,'second');
    
    % If in daylight savings time of that year, add an hour
    if is_DST(datestr(dateSerial,'mm/dd/yyyy')), dateSerial = dateSerial + datenum(0,0,0,1,0,0); end;
    
    % String version of date
    dateString = datestr(dateSerial); 
end