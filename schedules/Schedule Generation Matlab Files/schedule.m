% =========================================================================
% This script creates a fan schedule for REGCAP. The fan is told to turn on
% or off (1 or 0) for every minutes of a calendar year. Output is in the
% format of a 
% file
%
% Will Turner Mar 2011 LBNL
% =========================================================================
clear all
clc

% First day of year Jan 1st is a Sunday
%schedule = zeros(525600,1, 'uint32'); % Number of minutes in a year
schedule = zeros(525600,5, 'uint32'); % Number of minutes in a week

for i = 1:8760;
    for j = 1:60;
        schedule(4,i) = j;
    end
end
    

%weekOp = [1 0 0 1 0 0 0];

%dryerOnTime = 12.0;
%dryerOffTime = 15.0;

%day = 1;
%hour = 1;
%min = 1;
%day = 1;
%week = 1;

%for i = 1:10080;
%    if min == 60
%        hour = hour +1;
%        min = 1;
%    end
%    
%    if hour == 24
%        day = day +1;
%        hour = 1;
%    end
%    
%    if day == 7
%        week = week +1;
%        day = 1;
%    end
    
    %if weekOp(dayOfWeek)
    %    schedule(min) = 1;
    %end
%    min = min + 1;
%end
%minCount = 1;
%hourCount = 1;
%dayCount = 1;
