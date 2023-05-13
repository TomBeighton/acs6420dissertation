function calculateError(inputDen,actualDen,sensor1Den)
%calculateError.m created by Tom Beighton for ACS6420 Advaned Project
%   Function to take predicted values and generate RMSE and plots showing
%   different between Actual and Predicted Densities for sensor 2

%Plot densities
smoothingFactor = 20; %How much the data is smoothed to remove noise 
figure("Color","White","Name","Prediction Plotted")

%Conversion to hours
inputDen(:,1) = inputDen(:,1)/60; %Convert input data to hours
actualDen(:,1) = actualDen(:,1)/60; %Convert Actual data to hours
%Plot Results
plot(inputDen(:,1),movmean(inputDen(:,2),smoothingFactor),"Color","Red")
hold on
plot(actualDen(:,1),movmean(actualDen(:,2),smoothingFactor),"Color","Blue")
% plot(sensor1Den(:,1)/60,movmean(sensor1Den(:,2),smoothingFactor),"Color","Green")
% %(Above) Uncomment if you wish to plot the Sensor 1 Densities as well

%Format look of plot
xlabel("TimeSteps (hours)","FontSize",12,"FontWeight","bold","FontName","Times New Roman")
ylabel("Density (veh/km)","FontSize",12,"FontWeight","bold","FontName","Times New Roman")
legend("Predicted Sensor 2 Densities","Actual Sensor 2 Densities","FontSize",12,"FontName","Times New Roman")
h = gca;
h.LineWidth = 1.5;
xticks(0:2:24)

% %Calculate RMSE for each section of the day and overall:
TotalRMSE = sqrt(sum((actualDen(:,2) - inputDen(:,2)).^2)/length(inputDen(:,1)));
NightRMSEMorn = sqrt(sum((actualDen(actualDen(:,1)<6.5,2) - inputDen(inputDen(:,1)<6.5,2)).^2)/length(inputDen(inputDen(:,1)<6.5,1)));
NightRMSEEve = sqrt(sum((actualDen(actualDen(:,1)>18.5,2) - inputDen(inputDen(:,1)>18.5,2)).^2)/length(inputDen(inputDen(:,1)>18.5,1)));
MorningRushRMSE = sqrt(sum((actualDen((6.5<actualDen(:,1))&(actualDen(:,1)<10.5),2) - inputDen((6.5<actualDen(:,1))&(actualDen(:,1)<10.5),2)).^2)/length(inputDen((6.5<actualDen(:,1))&(actualDen(:,1)<10.5),2,1)));
MiddayRMSE = sqrt(sum((actualDen((10.5<actualDen(:,1))&(actualDen(:,1)<15.5),2) - inputDen((10.5<actualDen(:,1))&(actualDen(:,1)<15.5),2)).^2)/length(inputDen((10.5<actualDen(:,1))&(actualDen(:,1)<15.5),2,1)));
EveningRushRMSE = sqrt(sum((actualDen((15.5<actualDen(:,1))&(actualDen(:,1)<18.5),2) - inputDen((15.5<actualDen(:,1))&(actualDen(:,1)<18.5),2)).^2)/length(inputDen((15.5<actualDen(:,1))&(actualDen(:,1)<18.5),2,1)));
% %Display SSE Value to the Command line
message = ['Total RMSE: ',num2str(TotalRMSE), newline, 'Night Time Morning: ',num2str(NightRMSEMorn),newline,'Morning Rush: ', num2str(MorningRushRMSE),newline, 'Midday: ', num2str(MiddayRMSE),newline,'Evening Rush: ',num2str(EveningRushRMSE),newline,'Night Time Evening: ',num2str(NightRMSEEve)];
disp(message);

end

