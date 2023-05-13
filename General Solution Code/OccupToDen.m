function [sensor1Den,sensor2Den] = OccupToDen(sensor1Occup,sensor2Occup)
%OccupToDen.m Created by Tom Beighton for ACS6420 Advanced Project.
%   Function to convert Sensor 1 and Sensor 2 Occupancy Values into desnity
%   and extract timsteps and convert them into minutes.

%% Sensor 1
inputOccupancy = sensor1Occup;
inputhours = hours(inputOccupancy.TIME_HOURS);
mins = minutes(inputOccupancy.TIME_MINUTES);
timeHours = inputhours + mins;
timeMinutes = minutes(timeHours);
timeSeconds = seconds(timeHours);

%Occupancy conversion
occupancy = inputOccupancy.("occup.occupancy");
densities = ((1000/4.3))*(occupancy/100); %Density in VPKm

%Save output file
sensor1Den = [timeMinutes,densities];

%% Sensor 2
inputOccupancy = sensor2Occup;
inputhours = hours(inputOccupancy.TIME_HOURS);
mins = minutes(inputOccupancy.TIME_MINUTES);
timeHours = inputhours + mins;
timeMinutes = minutes(timeHours);
timeSeconds = seconds(timeHours);

%Occupancy conversion
occupancy = inputOccupancy.("occup.occupancy");
densities = (1000/4.3)*(occupancy/100); %Density in VPKm

%Save output file
sensor2Den = [timeMinutes,densities];

end

