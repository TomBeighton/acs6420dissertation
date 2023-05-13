%% Main Script to run multiple functions to process an extracted CSV and output error and simulation
% Created by Tom Beighton for ACS6420 Advanced Project
% To run please refer to README.Txt in the repository. Otherwise select
% run to output a simulation over 24 hours and calculate the error between
% predicted and actual values at Sensor 2. For a stretch of Penistone Road,
% Sheffield. 

close all
OutputToFDs = 0; %1 = Simulate , 0 = Plot FDs or Get FD parameters.
csvFilenameTest = "1DayPenistoneStretch1.csv"; %File to simulate
csvFilenameTrain = "3WeeksPenistoneStretch1.csv"; %File to tune FDs
sensor1name = "[SCC]1VHD1"; %Name of Sensor 1
plotSimulation = 1; %Where to plot the Time-Distance-Density graph
fundamentalDiagram = "Greenshields"; %What FD to use for simualtion or plotting
usage = "Parameters"; %"Parameters" = Get FD parameters, "Plot" = Plot FD
saveResult = 1; %Whether to save the graphs to the current workspace
trainingData = "Everyday"; %When tuning FDs what amount of data to use, can use "Mondays" or "Everyday"
step = 1; %There rarefaction wave step size.

if OutputToFDs == 1 %Check whether user wants to simulate or tune/plot FDs

    % Extract Flow and Occupancy Values
    [sensor1Flow,sensor1Occup,sensor2Flow,sensor2Occup] = ExtractCSV(csvFilenameTrain,sensor1name);
    [sensor1,sensor2] = MatchValues(sensor1Flow,sensor1Occup,sensor2Flow,sensor2Occup,trainingData);
    % (Optional) Output Data to FDs
    if usage == "plot"
        plotFDs(usage,fundamentalDiagram,sensor1,sensor2)
    else %Else determine FD parameters
        [rhoJam, QMax,rhoCrit,v0 ]  = plotFDs("parameters","CustomOverall",sensor1,sensor2);
    end

else
    % Extract Flow and Occupancy Values
    [sensor1Flow,sensor1Occup,sensor2Flow,sensor2Occup] = ExtractCSV(csvFilenameTest,sensor1name);

    %Convert Occupancy to Density and Extract Time Steps
    [sensor1Den,sensor2Den] = OccupToDen(sensor1Occup,sensor2Occup);

    %Sync up timesteps with simulation
    sensor1Den = [(sensor1Den(:,1) - sensor1Den(1,1)),sensor1Den(:,2)];
    sensor2Den = [(sensor2Den(:,1) - sensor2Den(1,1)),sensor2Den(:,2)];

    %INPUT TUNED FD PARAMETERS 
    %Penistone Road Stretch 1
    rhoJam = 200.5428; %veh/km
    QMax = 33.8031; %veh/min
    rhoCrit = 100.3718; %veh/km
    v0 = 0.3364; %km/min

%     %Penistone Road Stretch 2
%     rhoJam = 203.6461; %veh/km
%     QMax = 36.5997; %veh/min
%     rhoCrit = 101.9250; %veh/km
%     v0 = 0.3591; %km/min

    %Check whether an Empirical Fundamental Diagram is selected
    if fundamentalDiagram == "CustomOverall" || fundamentalDiagram == "CustomSensor1" || fundamentalDiagram == "CustomSensor2" 
        %Create Flow Density relationship fit and pass parameters into
        %simulation
        [sensor1Flow,sensor1Occup,sensor2Flow,sensor2Occup] = ExtractCSV(csvFilenameTrain,sensor1name);
        [sensor1,sensor2] = MatchValues(sensor1Flow,sensor1Occup,sensor2Flow,sensor2Occup,trainingData);
        [p1,p2] = plotFDs("coeff",fundamentalDiagram,sensor1,sensor2);
    else
        p1 = 0;
        p2 = 0;
    end
    
    %Run main simulation script to determine predicted values
    [predictedDen,bigZ] = simulation(sensor1Den,fundamentalDiagram,plotSimulation,rhoJam,QMax,rhoCrit,v0,step,p1,p2);

    % Error Calculation
    calculateError(predictedDen,sensor2Den,sensor1Den);
    if saveResult == 1 %Save error graph to current folder
        filename = sprintf('results%s.jpg',fundamentalDiagram);
        exportgraphics(gcf,filename)
    end
    
    %Save results to workspace if needed.
    save("OutputValues.mat", "predictedDen","sensor2Den");
end




