# acs6420dissertation
Code used in ACS6420 Final Year Project

==================================================================================================================================
Overview of Files:

This repository includes:

Folder: Data: #Folder containing all data needed for plotting of simulation
  Folder: Penistone Road Stretch 1
    1DayPenistoneStretch1.csv # 24 Hours of data from Penistone Road Stretch 1
    3WeeksPenistoneStretch2.csv # 3 Weeks of data from Penistone Road Stretch 1
  Folder: Penistone Road Stretch 2 
    1DayPenistoneStretch1.csv # 24 Hours of data from Penistone Road Stretch 1
    3WeeksPenistoneStretch2.csv # 3 Weeks of data from Penistone Road Stretch 1
Folder: General Solution Code
    simulation.m #Function containing the LWR Model General Solution and outputs Time-Distance-Density graph
    ExtractCSV.m #Function to extract input data from csv into Flow and Occupancy Data
    MatchValues.m #Function to convert Occupancy into Density and match Flow measurements with density measurements for Both Sensors
    OccupToDen #Function to convert Occupancy into Density and extract timsteps of recording for both sensors.
    calculateError.m #Function to calculate RMSE for all sections of the day and plot graph of predicted and actual Sensor 2 values
    linexline.m #External function to calculate the intercept of each shockwave.
    FDs.m #Function to calculate the value of Lambda based upon the Fundamental Diagram specified by the user
    plotFDs.m #Function to plot Empirical Fundamental Diagrams 
    ProcessSimulation.m #Main .m file whihc when run calls previously described functions to achieve simulation or plot specified.

====================================================================================================================================
Operation:

Example: Plot the simulation of a stretch of road 100m long using 1 Day of Sensor Data and Tuned Greenshields FD using 3 weeks of data.

Step 1: Open MATLAB and open ProcessSimulation.m and Simulation.m ensure that the data for simulation and for tuning is added to MATLAB's path so it can be accessed
Step 2: Prepare the simulaion to run. In ProcessSimulation.m enter the correct file names for testing and tuning data, the correct name of sensor 1 (Found in csv), and specify the step size.
Step 3: Other options include: Training with different amounts of data "Mondays" or "Weekdays"
Step 4: Now tune the Fundamental Diagram by ensuring the fundamentDiagram == "CustomOverall" and usage=="Parameters" and most importantly, OuputToFDs == 1
Step 5: Run the program. The workspace will now be populated with variables [rhojam,rhocrit,v0,Qmax]
Step 6: Enter tuned variables on appropriate lines in ProcessSimulation.m and now change fundamentalDiagram == "Greenshields", OuputsToFDs == 0 and plotSimulation == 1
Step 7: If you wish to save the graph of the results automatically please set saveResult == 1;
Step 8: Open Simulation.m and adjust the endDistance to the distance between sensors on your target stretch and vf to the value found from literature. Note velocity is in km/min
Step 9: Open ProcessSimulation.m and you can now run the program. This might take some time. In the command window RMSE values will be shown and graphs displayed for Predicted vs Actual and Time-Distance-Density.

You can now run the simulation again with a different Fundamental Diagram by changing the fundamentalDiagram variable to any of the following options:
["CustomOverall","CustomSensor1","CustomSensor2","Greenshields","TrafficFreeFlow","Underwood","Greenberg","Newell-Daganzo"]

NOTE: The current as deposited in this repository is able to perform the above example simply by selected run provided the data is added to the path. No values need to be changed but above gives some steps if you wish to use your own.

======================================================================================================================================
Data:

Data can be extracted from the Sheffield Urban Flows Observatory using the following link on a University LAN connected desktop:
https://ufdev21.shef.ac.uk/sufobin/sufoDXT

You can extract data from any TWO sensors that have both FLOW and OCCUPANCY measurements (Make sure to extract ONLY the Flow and Occupancy). Make sure to select the amount of time and do not select too much data to extract otherwise extraction and simulation will take a long time and is not guarenteed to work.

NOTE: that the data is in a minimal format suitable for MATLAB (one of the dropdown options).

Thanks for using my code!
