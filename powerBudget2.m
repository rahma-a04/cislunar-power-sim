function powerBudget()
%This function calculates the power drain and compares it to the power
%generated in order to determine the state of charge of the battery over
%the mission lifetime
clear all
close all
clc

%Set defaults for all the graphs
set(0,'DefaultAxesFontSize',18,'DefaultTextFontSize',18)

%Load the mission data
S = load('SLSLunarTraj_mar25.mat');
%Define the useful parameters from this data
%tarray = S.tarray
tarray = S.tarray2;
indicearray = S.indicearray2;
t = S.timespans2;

%Find the time stamps of the different modes
stages = zeros(length(t),1);
stages(1) = t(1);
for ii = 2:length(t)
    stages(ii) = stages(ii-1) + t(ii);
end

%The time between each maneuver
dtPulse = [54000 56970 64470 74110];

%Number of divisions of data
n = length(tarray);

%Find the change in time between steps
for i = 1:n-1
    dt(i) = tarray(i+1)-tarray(i);
end

dt = [dt dt(end)];

%Define Power Modes
postSep = 1;
init = 2;
downlink = 3;
reorient = 4;
elect = 5;
cruise = 6;
fire = 7;
seq = 8; %Fire, electrolyze, reorient sequence

%% THIS SECTION NEEDS TO BE CHANGED %%

%This is the expected maximum duration of each mode in seconds.  These 
%values should be updated as new estimates become available
modeDur = [60 2760 60 300 5097 1 3];

%This vector is the expected power draw at each mode.  These values should
%be updated as the power budget gets finalized.
powerDraw(1:7) = [4.77 
                    0.26 
                    9.19 
                    4.34 
                    8.15
                    3.48 
                    1.68]';
       
powerDraw(8) = (powerDraw(elect)*modeDur(elect) + ...
    powerDraw(reorient)*modeDur(reorient) + ...
    powerDraw(fire)*modeDur(fire))/...
    (modeDur(elect) + modeDur(fire) + modeDur(reorient));
powerDraw(9) = (7.38*5697 + ...
    powerDraw(reorient)*modeDur(reorient) + ...
    powerDraw(fire)*modeDur(fire))/...
    (5697 + modeDur(fire) + modeDur(reorient));
powerDraw(10) = (6.99*6049 + ...
    powerDraw(reorient)*modeDur(reorient) + ...
    powerDraw(fire)*modeDur(fire))/...
    (6049 + modeDur(fire) + modeDur(reorient));
powerDraw(11) = (5.83*7411 + ...
    powerDraw(reorient)*modeDur(reorient) + ...
    powerDraw(fire)*modeDur(fire))/...
    (7411 + modeDur(fire) + modeDur(reorient));

%% Find the power usage at each time step from the power modes

%Initialize the vector that defines the power mode at each time step
powerMode = zeros(n,1);

%Determined the power mode at each time step
for i=1:n
    if tarray(i) <= modeDur(1)
        powerMode(i) = postSep;
    elseif tarray(i) <= (modeDur(1) + modeDur(2))
        powerMode(i) = init;
    elseif tarray(i) <= (modeDur(1) + modeDur(2) + modeDur(3))
        powerMode(i) = downlink;
    elseif tarray(i) <= (modeDur(1) + modeDur(2) + modeDur(3) + modeDur(4))
        powerMode(i) = reorient;
    elseif ~isempty(find(tarray(indicearray) == tarray(i)))
        if tarray(i) < 2e6
            powerMode(i) = 10;
            time = tarray(i);
            k = 1;
            while tarray(i+k) <= time+dtPulse(3)
                powerMode(i+k) = 10;
                k=k+1;
            end
        elseif tarray(i) < 4e6
            powerMode(i) = 10;
            time = tarray(i);
            k = 1;
            while tarray(i+k) <= time+dtPulse(3)
                powerMode(i+k) = 10;
                k=k+1;
            end
        elseif tarray(i) < 6e6 %CHANGED
            powerMode(i) = 11;
            time = tarray(i);
            k = 1;
            %disp(time+dtPulse(4))
            while tarray(i+k) <= time+dtPulse(4)
                powerMode(i+k) = 11;
                k=k+1;
            end
        end
    end
end

%If the satellite isn't defined in one of the modes above, it is cruising
powerMode(powerMode==0) = 6;

%Define the associated power draw expected at each time step now that the
%modes have been determined for each time step
for i = 1:n
    powerUsage(i) = powerDraw(powerMode(i));
end

%Random vector to plot the stages on the graphs
Y = zeros(length(stages));

% Plot the power used over the mission lifetime, and indicate the different
% stages of the satellite mission
figure
plot(tarray,powerUsage,stages,Y)
title('Power Usage over the Mission Life')
xlabel('Time (s)')
ylabel('Power (W)')
axis([0 tarray(end) 0 10])

%Call the powerCalc function to determine the power generated by the solar
%panels at each time step
powerIn = powerCalc2;

%The net power at each step is the power generated minus the power used
netPower = powerIn' - powerUsage;

%% THIS SECTION TO BE CHANGED %%

%The minimum battery capacity.  If the battery design changes, this value
%needs to be updated.
minBatCap = 207792; %[Ws]

%% Battery State Calculations
%Initalize the state of the battery charge vector
batCharge = zeros(n,1);

%Start with the battery fully charged
batCharge(1) = minBatCap;

%The battery charge at every step is the minimum between the previous state
%of the battery plus the net power over the previous time step, or the
%minimum battery capacity which is the battery fully charged
for i = 2:n
    batCharge(i) = min(minBatCap,(batCharge(i-1) + netPower(i)*dt(i)));
end

%The percentage of the battery that is charged
batCharge = batCharge/minBatCap*100;

%Plot the net power over the mission lifetime
figure
plot(tarray,netPower,stages,Y,'x-','MarkerSize',10)
title('Net Power over the Mission Life')
xlabel('Time (s)')
ylabel('Power (W)')
axis([0 tarray(end) -2 4])
whitebg('k')

%Plot the state of the battery charge over the mission lifetime
figure
plot(tarray,batCharge)
title('State of Battery Charge over the Mission Life')
xlabel('Time (s)')
ylabel('Percent Charge (%)')
axis([0 tarray(end) 0 100])
whitebg('k')


end