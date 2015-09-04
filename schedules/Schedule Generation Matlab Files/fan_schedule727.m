% =========================================================================
% This script creates a fan schedule for REGCAP. The fan is told to turn on
% or off (1 or 0) for every minute of a calendar year by calling the
% 'fanSwitch' function. Output is in the 'schedule' array.
%
% Will Turner Mar 2011 LBNL
% =========================================================================

% First day of year Jan 1st is a Sunday

schedule = zeros(525600,8, 'uint8'); % Number of minutes in a week
%schedule = zeros(525600,8); % Number of minutes in a week

day = 1;    % schedule column 1
hour = 0;   % schedule column 2
min = 0;    % schedule column 3
% dryer       schedule column 4
% kitchen fan schedule column 5
% bath fan    schedule column 6


% Create a year
for i = 1:525600;
    if min == 60
        min = 0;
        hour = hour +1;
    end
    
    if hour == 24
        hour = 0;
        day = day +1;
    end
    
    if day == 8
        day = 1;
    end
    
    schedule(i,1) = day;
    schedule(i,2) = hour;
    schedule(i,3) = min;
    
    min = min +1;
end

%         [S M T W T F S]   1 = ON, 0 = OFF
week    = [0 1 1 1 1 1 0];
weekend = [1 0 0 0 0 0 1];

%% =========================== DRYER =======================================
nFan = 4;
dryer.days = [1 0 0 1 0 0 0];
dryer.on = [11 00];     % Dryer start time in 24hr notation
dryer.mins = 180;       % Number of consecutive hours of dryer operation

schedule = fanSwitch(schedule, nFan, dryer.days, dryer.on, dryer.mins);

%% ========================= KITCHEN =======================================

% AM weekends
nFan = 5;
kitchenAM.on = [9 30];
kitchenAM.mins = 30;

schedule = fanSwitch(schedule, nFan, weekend, kitchenAM.on, kitchenAM.mins);

% PM every day
nFan = 5;
kitchenPM.on = [19 30];
kitchenPM.mins = 60;
 
schedule = fanSwitch(schedule, nFan, (week + weekend), kitchenPM.on, kitchenPM.mins);

%% ========================= BATHROOM ======================================

% BATHROOM 1: WEEK DAYS
occupants = 2;
nFan = 6;
on = [6 0];
mins = 60;

wcStart = 19;
wcSpan = 5;
wcMins = 10;

for a = 1:365;
    schedule = fanWeek(schedule, nFan, week, on, mins, wcStart, wcSpan, wcMins, occupants, a);
end

% BATHROOM 1: WEEKENDS
showerStart = 7;
showerSpan = 12;
showerMins = 30;

wcStart = 7;
wcSpan = 16;
wcMins = 10;
for a = 1:365;
    schedule = fanWeekend(schedule, nFan, weekend, showerStart, showerSpan, showerMins, wcStart, wcSpan, wcMins, occupants, a);
end
    % -------------------------------------------------------------------------
%% BATHROOM 2: WEEK DAYS
occupants = 2;
nFan = 7;
on = [6 0];
mins = 60;

wcStart = 19;
wcSpan = 5;
wcMins = 10;

for a = 1:365;
    schedule = fanWeek(schedule, nFan, week, on, mins, wcStart, wcSpan, wcMins, occupants, a);
end
% BATHROOM 2: WEEKENDS
showerStart = 7;
showerSpan = 12;
showerMins = 30;

wcStart = 7;
wcSpan = 16;
wcMins = 10;

for a = 1:365;
    schedule = fanWeekend(schedule, nFan, weekend, showerStart, showerSpan, showerMins, wcStart, wcSpan, wcMins, occupants, a);
end
% -------------------------------------------------------------------------
%% BATHROOM 3: WEEK DAYS
occupants = 1;
nFan = 8;
on = [6 30];
mins = 30;

wcStart = 19;
wcSpan = 5;
wcMins = 10;

for a = 1:365;
    schedule = fanWeek(schedule, nFan, week, on, mins, wcStart, wcSpan, wcMins, occupants, a);
end

% BATHROOM 3: WEEKENDS
showerStart = 7;
showerSpan = 12;
showerMins = 30;

wcStart = 7;
wcSpan = 16;
wcMins = 10;

for a = 1:365;
    schedule = fanWeekend(schedule, nFan, weekend, showerStart, showerSpan, showerMins, wcStart, wcSpan, wcMins, occupants, a);
end