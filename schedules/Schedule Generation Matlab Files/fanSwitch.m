function schedule = fanSwitch(schedule, nFan, days, startTime, minutes)
% FAN SWITCH function
% Switches fans on (1) or off (0) for one calendar year, minute-by-minute

for i = 1:525600;
    if days(schedule(i,1))
        if schedule(i,2) == startTime(1) && schedule(i,3) == startTime(2)
            for j = 0:minutes-1
                schedule(i+j,nFan) = 1;
            end
            i = i + minutes;
        end
    end
end