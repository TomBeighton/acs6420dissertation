function [sensor1FlowValues,sensor1OccupValues,sensor2FlowValues,sensor2OccupValues] = ExtractCSV(CSVfilename,sensor1Name)
%ExtractCSV.m Created by Tom Beighton for ACS6420 Advanced Project
%   Function to take input data and extract Flow and Occupancy values for
%   Sensor 1 and Sensor 2

inputCSV = readcell(CSVfilename); %Read in filename

%Find Flow Headers
indexFlowHeaders = find(strcmp(inputCSV(:,1),'flows.time'));

%Find Occupancy Headers
indexOccupHeaders = find(strcmp(inputCSV(:,1),'occup.time'));

%ExtractValues into Table
flowValues1 = cell2table(inputCSV(indexFlowHeaders(1) + 1:indexFlowHeaders(2) - 3,:),'VariableNames',inputCSV(indexFlowHeaders(1),:));
flowValues2 = cell2table(inputCSV(indexFlowHeaders(2) + 1:indexOccupHeaders(1) - 3,:),'VariableNames',inputCSV(indexFlowHeaders(2),:));
occupValues1 = cell2table(inputCSV(indexOccupHeaders(1) + 1:indexOccupHeaders(2) - 3,:),'VariableNames',inputCSV(indexOccupHeaders(1),:));
occupValues2 = cell2table(inputCSV(indexOccupHeaders(2) + 1:length(inputCSV) - 1,:),'VariableNames',inputCSV(indexOccupHeaders(2),:));



%Split into Sensor 1 and Sensor 2 Values
if flowValues1.("flows.sensor")(1) == sensor1Name
    sensor1FlowValues = flowValues1;
    sensor2FlowValues = flowValues2;
else
    sensor1FlowValues = flowValues2;
    sensor2FlowValues = flowValues1;
end
if occupValues1.("occup.sensor")(1) == sensor1Name
    sensor1OccupValues = occupValues1;
    sensor2OccupValues = occupValues2;
else
    sensor1OccupValues = occupValues2;
    sensor2OccupValues = occupValues1;
end



end

