function schedule = fanWeek(schedule, nFan, days, startTime, minutes, wcStart, wcSpan, wcMins, occupants, a)
% FAN SWITCH function
% Switches fans on (1) or off (0) for one calendar year, minute-by-minute
%% Set random start times
wcH = [(wcStart + randi(wcSpan,1)-1) (wcStart + randi(wcSpan,1)-1)];
wcM = [((randi(2)-1)*30) ((randi(2)-1)*30)];

if wcH(1) == wcH(2) && wcM(1) == wcM(2)
    while wcM(1) == wcM(2)
        wcM(2) = (randi(2)-1)*30;
    end
end
%% Morning regular schedule
for i = (a*1440-1439):(a*1440);
    if days(schedule(i,1))
        if schedule(i,2) == startTime(1) && schedule(i,3) == startTime(2)
            for j = 0:minutes-1
                schedule(i+j,nFan) = 1;
            end
            i = i + minutes;
        end
    end
end    
%% Evening variable schedule #1
for i = (a*1440-1439):(a*1440);
    if days(schedule(i,1))
        if schedule(i,2) == wcH(1) && schedule(i,3) == wcM(1)
            for j = 0:wcMins-1
                schedule(i+j,nFan) = 1;
            end
            i = i + wcMins;
        end
    end
end
%% Evening variable schedule #2
if occupants == 2;
    for i = (a*1440-1439):(a*1440);
        if days(schedule(i,1))
            if schedule(i,2) == wcH(2) && schedule(i,3) == wcM(2)
                for j = 0:wcMins-1
                    schedule(i+j,nFan) = 1;
                end
                i = i + wcMins;
            end
        end
    end
end