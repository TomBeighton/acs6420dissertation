function [Sensor1Matching,Sensor2Matching] = MatchValues(sensor1Flow,sensor1Occup,sensor2Flow,sensor2Occup,inclusions)
%MatchValues.m Created by Tom Beighton for ACS6420 Advanced Project.
%       Function to take Sensor 1 and Sensor 2 values, convert Occupancy
%       Values to density, and match Flow and Density values recorded at
%       the same time. 
   
    %Sensor 1
    Flow = sensor1Flow;
    Occupancy = sensor1Occup;


    if inclusions == "Mondays" %Extract only Mondays 
        Flow = Flow(strcmp(Flow.TIME_DOW,'Mo'),:);
        Occupancy = Occupancy(strcmp(Occupancy.TIME_DOW,'Mo'),:);
    elseif inclusions == "Weekdays" %Extract only weekdays
        weekdays = {'Mo','Tu','We','Th','Fr'};
        Flow = Flow((not(strcmp(Flow.TIME_DOW,'Sa')))&(not(strcmp(Flow.TIME_DOW,'Su'))),:);
        Occupancy = Occupancy((not(strcmp(Occupancy.TIME_DOW,'Sa')))&(not(strcmp(Occupancy.TIME_DOW,'Su'))),:);
    end

    FlowMat = [Flow.TIME_MINUTES,Flow.TIME_HOURS,Flow.TIME_DAY,Flow.("flows.flow")];
   
    OccupMat = [Occupancy.TIME_MINUTES,Occupancy.TIME_HOURS,Occupancy.TIME_DAY,Occupancy.("occup.occupancy")];
    %Ensure unique values of all Flow and Occupancy Recordings:
    [FlowMatDates,ia,~] = unique(FlowMat(:,1:end-1),"rows",'stable');
    FlowMat = FlowMat(ia,:);
    [~,ia,~] = unique(OccupMat(:,1:end-1),"rows",'stable');
    OccupMat = OccupMat(ia,:);
    
 

    %Find matching records for both datasets by searching through each
    MatchingOccup = OccupMat(ismember(OccupMat(:,1:end-1),FlowMatDates,"rows"),end);
    MatchingFlow = FlowMat(ismember(FlowMat(:,1:end-1),OccupMat(:,1:end-1),"rows"),end);
    matchingTimes = [MatchingOccup, MatchingFlow]; %Join together matching times into one matrix.
    

    occupancies = matchingTimes(:,1);
    y = matchingTimes(:,2);

    %Convert Occupancies to Densities
    x = ((1000/4.3)/100)*occupancies; %Density in VPKm

    %SaveValues
    Sensor1Matching = [x,y];

    %Sensor 2
    Flow = sensor2Flow;
    Occupancy = sensor2Occup;


     if inclusions == "Mondays" %Extract only Mondays 
        Flow = Flow(strcmp(Flow.TIME_DOW,'Mo'),:);
        Occupancy = Occupancy(strcmp(Occupancy.TIME_DOW,'Mo'),:);
    elseif inclusions == "Weekdays" %Extract only weekdays
        weekdays = {'Mo','Tu','We','Th','Fr'};
        Flow = Flow((not(strcmp(Flow.TIME_DOW,'Sa')))&(not(strcmp(Flow.TIME_DOW,'Su'))),:);
        Occupancy = Occupancy((not(strcmp(Occupancy.TIME_DOW,'Sa')))&(not(strcmp(Occupancy.TIME_DOW,'Su'))),:);
    end

    FlowMat = [Flow.TIME_MINUTES,Flow.TIME_HOURS,Flow.TIME_DAY,Flow.("flows.flow")];
   
    OccupMat = [Occupancy.TIME_MINUTES,Occupancy.TIME_HOURS,Occupancy.TIME_DAY,Occupancy.("occup.occupancy")];


    %Ensure unique values of all Flow and Occupancy Recordings:
    [FlowMatDates,ia,~] = unique(FlowMat(:,1:end-1),"rows",'stable');
    FlowMat = FlowMat(ia,:);
    [~,ia,~] = unique(OccupMat(:,1:end-1),"rows",'stable');
    OccupMat = OccupMat(ia,:);


    %Find matching records for both datasets by searching through each
    MatchingOccup = OccupMat(ismember(OccupMat(:,1:end-1),FlowMatDates,"rows"),end);
    MatchingFlow = FlowMat(ismember(FlowMat(:,1:end-1),OccupMat(:,1:end-1),"rows"),end);
    matchingTimes = [MatchingOccup, MatchingFlow]; %Join together matching times into one matrix.

    %Convert Occupancy to Density.
    occupancies = matchingTimes(:,1);
    y = matchingTimes(:,2);
    x = ((1000/4.3)/100)*occupancies; %Density in VPKm

    %SaveValues
    Sensor2Matching = [x,y];


end

