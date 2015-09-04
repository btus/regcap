function schedule = fanWeekend(schedule, nFan, days, showerStart, showerSpan, showerMins, wcStart, wcSpan, wcMins, occupants, a)
% FAN SWITCH function
% Switches fans on (1) or off (0) for one calendar year, minute-by-minute
%%
showerH = [(showerStart + randi(showerSpan,1)-1) (showerStart + randi(showerSpan,1)-1)];  % Random integer between 7 and 18
showerM = [((randi(2)-1)*30) ((randi(2)-1)*30)];                      % Random integer 0 or 30
if showerH(1) == showerH(2) && showerM(1) == showerM(2)
    while showerM(1) == showerM(2)
        showerM(2) = (randi(2)-1)*30;
    end
end
%%
wcH = [(wcStart + randi(wcSpan,1)-1) (wcStart + randi(wcSpan,1)-1)];
if wcH(1) == showerH(1) || wcH(1) == showerH(2)
    while wcH(1) == showerH(1) || wcH(1) == showerH(2)
        wcH(1) = wcStart + randi(wcSpan,1)-1;
    end
end
%%
if wcH(2) == showerH(1) || wcH(2) == showerH(2) || wcH(1) == wcH(2)
    while wcH(2) == showerH(1) || wcH(2) == showerH(2) || wcH(1) == wcH(2)
        wcH(2) = wcStart + randi(wcSpan,1)-1;
    end
end

wcM = [((randi(2)-1)*30) ((randi(2)-1)*30)];
%%
for i = (a*1440-1439):(a*1440);
    if days(schedule(i,1))
        if schedule(i,2) == showerH(1) && schedule(i,3) == showerM(1)
            for j = 0:showerMins-1
                schedule(i+j,nFan) = 1;
            end
            i = i + showerMins;
        end
    end
end
%%
if occupants == 2;
    for i = (a*1440-1439):(a*1440);
        if days(schedule(i,1))
            if schedule(i,2) == showerH(2) && schedule(i,3) == showerM(2)
                for j = 0:showerMins-1
                    schedule(i+j,nFan) = 1;
                end
                i = i + showerMins;
            end
        end
    end
end
%%
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
%%
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