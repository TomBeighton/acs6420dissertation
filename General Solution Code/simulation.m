function [outputDen,bigZ] = simulation(inputData,fundamentalDiagram,plotdecision,rhoMax,QMax,rhoCrit,v0,step,p1,p2)
%simulation.m Created by Tom Beighton for ACS6420 Advanced Project
%Takes input parameters from ProcessSimulation.m and uses General Solution
%of LWR model to plot a Time-Distance-Density graph from which predicted
%sensor 2 values are determined.

%% Initilise Simulation
%Define variables for the Fundamental Diagram to be used:
vf = 0.672; %Traffic Free Flow Speed km/min
FD = fundamentalDiagram; %Choose the Fundamental Diagram

%Define Input Densities:
initialdensity = inputData(1,2); %Extract Sensor 1 data from input
inputDen = inputData(2:end,2);
timeSteps = inputData(2:end,1);
resolution = 0.01; %Choose resolution of graph, Warning high res= long running times.
endDistance = 0.1; %Kms 
% endDistance = 2;
%Define piecewise constants between Sensor 1 and Sensor 2
den = initialdensity; %Densities between sensor 1 and 2 for initlisation.
bounds = endDistance; %Upperbounds of density regions

%Start program
carryIntercept = 0; %Intercept that makes bottom of meshgrid. Starts at the first timestep
endTime = timeSteps(end) + (timeSteps(2)-timeSteps(1)); %Sets End time to one step above last entry.
nextDenTime = timeSteps(1); %The next timestep
nextDenIndex = 2; %Index of the next timestep within the inputDen array
calculateIntercepts= 1;


bigZ = zeros(1/resolution*(length(timeSteps)+1),endDistance/resolution); %Large matrix of zeros to hold overall graph
[bigX,bigT] = meshgrid(linspace(0,endDistance,endDistance/resolution),linspace(inputData(1,1),endTime,height(bigZ)));
% Main Loop
while carryIntercept < endTime %Iterate until end of simulation time.

    %Check if current time is equal to time to insert next timestep
    if nextDenIndex > length(inputDen)
        %Carry on simulation
        nextDenTime = endTime;
        calculateIntercepts = 1;
    end
    if ismember(carryIntercept,timeSteps) == 1 %Check if next time step has been reached
        calculateIntercepts = 1;
        if nextDenTime ~= endTime %Check if new timestep is not the last one
            nextDenTime = timeSteps(nextDenIndex); %Update next entry time
        end
        nextden = inputDen(nextDenIndex - 1); %Get the next density
        nextDenIndex = nextDenIndex + 1;

        %Appends any new densities to the front of the density list
        if nextden < den(1) %Check for each discontinuity condition.
            %Shockwave
            den = [nextden,den];
            bounds = [0,bounds]; %Append density and boundary to front of den matrix
        elseif nextden > den(1)
            %Rarefaction
            bounds = [zeros(1,length(nextden:-step:den(1))-1),bounds];
            den = [nextden:-step:den(1)+step,den];
        else
            %Contact Discontinuity nothing happens
        end
    end
    if length(den) == 1 %Check if the first density is selected
        [t,x] = meshgrid(linspace(0,endDistance,endDistance/resolution),linspace(0,nextDenTime-carryIntercept,(nextDenTime-carryIntercept)/resolution));
        z = ones(height(t),width(t)) * den(1);
        xIndexes = (carryIntercept/resolution +1):1:(carryIntercept/resolution + height(z) );
        carryIntercept = nextDenTime;
        calculateIntercepts = 0;
    end

    if calculateIntercepts ==1 %Calculate intercepts of shockwaves
        %Generate line gradients (Lambdas)
        Lambdas = [;];
        x=linspace(-100,endDistance+100,1000);
        for i = 2: length(den)

            rhoPlus = den(i);
            rhoMinus = den(i-1);
            %Runs function to calculate Lambda based on input Fundamental
            %Diagram
            lambda = Fds(rhoMinus,rhoPlus,FD,vf,v0,rhoMax,QMax, rhoCrit,p1,p2);
            lambda = round(lambda,4);
            if sum(isnan(lambda)) > 0 %Checks for contact discontinuity
                lambda = 0;
            end

            try
             Lambdas = [Lambdas,[lambda;bounds(i-1);den(i);den(i-1)]]; 
             %Creates an array of Lambda values, their intercept with x-axis and densities either side of shockwave.
            catch
                
            end
        end
        if width(Lambdas) == 1
            %Only two densities and therefore no intersections possible
            %Now, plot Lambda to the next timestep

            minIntersect = [];

        else
            % Generates Intercepts of Lines
            intercepts = [;];
            for i = 1: width(Lambdas) -1
                x=linspace(-100,endDistance+100,1000);
                %Find intersect with Neighbour

                if Lambdas(1,i) == 0 %Calculates equation of shockwave line
                    x1 = x;
                    t1 = (x- Lambdas(2, i+1))/Lambdas(1,i+1);
                    t2 = x;
                    x2 = Lambdas(2,i)*ones(1,length(x));
                elseif Lambdas(1,i+1) == 0
                    t2 = (x- Lambdas(2, i))/Lambdas(1,i);
                    t1 = x;
                    x1 = Lambdas(2,i+1)*ones(2,length(x));
                    x2 = x;
                else
                    t1 = (x- Lambdas(2, i+1))/Lambdas(1,i+1);
                    t2 = (x- Lambdas(2, i))/Lambdas(1,i);
                    x1 = x;
                    x2 = x;
                end

                L1x = [x1(1),x1(end)];
                L1y = [t1(1),t1(end)];
                L2x = [x2(1),x2(end)];
                L2y = [t2(1),t2(end)];
                [xi,yi] = linexline(L1x,L1y,L2x,L2y,0); %Run external script to find intercept point.
                intercept = [xi;yi];
                if isempty(intercept) == 1
                    intercept(1,1) = 0;
                    intercept(2,1) = 0;
                end


                %Construct Intercept Matrix
                intercepts = [intercepts,[round(intercept,5);Lambdas(1,i+1);Lambdas(1,i);den(i);den(i+1);den(i+2)]];
            end

            %Check for any valid intercepts below the next input density.
            viableIntercepts = intercepts(2,:);
            minIntersect = min(viableIntercepts(viableIntercepts>resolution));
        end


        if isempty(minIntersect) == 1 ||  (minIntersect + carryIntercept > nextDenTime)%No intercepts above y = 0

            %Continue on lines till next timestep.

            minIntersect = nextDenTime-carryIntercept;
            % Create new bounds Matrix
            newZInterceptsX = [];
            for i =1 :width(Lambdas) %Find points of lines at intersect point
                if Lambdas(1,i) == 0 %If straight line next point is directly upwards
                    point = Lambdas(2,i);
                else
                    point = (minIntersect * Lambdas(1,i)+Lambdas(2,i));
                end
                if point > 0 && point < endDistance %Check whether intercept is within bounds of road stretch
                    newZInterceptsX = [newZInterceptsX,[point;Lambdas(4,i);Lambdas(3,i);1]];%X intercepts with horizontal line
                else
                    newZInterceptsX = [newZInterceptsX,[point;Lambdas(4,i);Lambdas(3,i);0]];%X intercepts with horizontal line
                end
            end

            if sum(newZInterceptsX(4,:))~=0 
                validDensities = newZInterceptsX(1:3,newZInterceptsX(4,:) == 1);
                newden = [validDensities(2,:),validDensities(3,end)]; %Define new densities and intercepts.
                bounds = [validDensities(1,:), endDistance];

            else
                %There are no points which intersect next time step between
                %0 - endDistance.
                distances = [];
                halfway = endDistance/2;
                for i = 1:width(newZInterceptsX) %Find the density that will dominate the graph.
                    if newZInterceptsX(1,i) <= 0
                        distances = [distances,[abs(newZInterceptsX(1,i)-halfway);newZInterceptsX(3,i)]];
                    else
                        distances = [distances,[newZInterceptsX(1,i)-halfway;newZInterceptsX(2,i)]];
                    end
                end
                [uniqueDist,ia,~] = unique(round(distances(1,distances(1,:)>0),4),'first');
                if length(uniqueDist) == length(distances(1,distances(1,:)>0))
                    newden = distances(2,distances(1,:) == min(distances(1,distances(1,:)>0)));
                else
                    newden = distances(2,distances(1,:) == min(distances(1,ia)));
                end
                if  width(newden) > 1
                    newden = newden(1);
                end

                bounds = endDistance;
            end

            %Remaining info is bounds continaing sections at next time step
            [t,x] = meshgrid(linspace(0,endDistance,endDistance/resolution),linspace(0,nextDenTime-carryIntercept,(nextDenTime-carryIntercept)/resolution ));

            if height(t) ~= (nextDenTime-carryIntercept*(1/resolution))
                [t,x] = meshgrid(linspace(0,endDistance,endDistance/resolution),linspace(0,nextDenTime-carryIntercept,(nextDenTime-carryIntercept)/resolution+1));
            end
            z = zeros(height(t),width(t));
            xIndexes = (round(carryIntercept/resolution) +1):1:(round(carryIntercept/resolution)+ height(z));
            carryIntercept = nextDenTime;



        else %Intercept is before next timestep

            index = find(viableIntercepts==minIntersect); %Extract intercept
            % Create new bounds Matrix
            newZInterceptsX = [];
            for i =1 :width(Lambdas) %Find points of lines at intersect point
                if Lambdas(1,i) == 0 %If straight line next point is directly upwards
                    point = Lambdas(2,i);
                else
                    point = (minIntersect * Lambdas(1,i)+Lambdas(2,i));
                end
                if point > 0 && point < endDistance
                    newZInterceptsX = [newZInterceptsX,[point;Lambdas(4,i);Lambdas(3,i);1]];%X intercepts with horizontal line
                else
                    newZInterceptsX = [newZInterceptsX,[point;Lambdas(4,i);Lambdas(3,i);0]];%X intercepts with horizontal line
                end
            end
            tempBounds = [];
            newden = [];

            if sum(newZInterceptsX(4,:))~=0
                %NewDen equal to all densities
                newden = [newZInterceptsX(2,:),newZInterceptsX(3,end)];
                tempBounds =  newZInterceptsX(1,:);

                if index == 1
                    newden = [newden(1),newden(index+2:end)];
                else
                    newden = [newden(1:index),newden(index+2:end)];
                end
                bounds = [unique(round(tempBounds,4),'stable'),endDistance];
            else
                %Intercept is located outside the range of 0-1 and
                %therefore only 1 density is carried forward.
                distances = [];
                halfway = endDistance/2;
                for i = 1:width(newZInterceptsX)
                    if newZInterceptsX(1,i) < 0
                        distances = [distances,[abs(newZInterceptsX(1,i)-halfway);newZInterceptsX(3,i)]];
                    else
                        distances = [distances,[newZInterceptsX(1,i)-halfway;newZInterceptsX(2,i)]];
                    end
                end
                [uniqueDist,ia,~] = unique(round(distances(1,distances(1,:)>0),4),'first');
                if length(uniqueDist) == length(distances(1,distances(1,:)>0))
                    newden = distances(2,distances(1,:) == min(distances(1,distances(1,:)>0)));
                else
                    newden = distances(2,distances(1,:) == min(distances(1,ia)));
                end
                if  width(newden) > 1
                    newden = newden(1);
                end

                bounds = endDistance;
            end
            % Draw Meshgrid of Densities up to intercept point
            [t,x] = meshgrid(linspace(0,endDistance,endDistance/resolution),linspace(0,intercepts(2,index),(ceil(intercepts(2,index)*100)/100)/resolution));
            z = zeros(height(t),width(t));
            xIndexes = (round(carryIntercept/resolution) +1):1:(round(carryIntercept/resolution) + height(z) );
            carryIntercept = carryIntercept + ceil(intercepts(2,index)*100)/100; %update current level
        end

        %Fill in meshgrid to create graph until next timestep or intercept
        %point
        for i = 2:width(Lambdas)
            if i == width(Lambdas)
                if  Lambdas(1,i) > 0 && Lambdas(1,i-1) < 0
                    indexes =  (x >= ((t- Lambdas(2, i))/Lambdas(1,i)) & x>= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                    z(indexes) = den(i);
                    indexes = (x <= (t- Lambdas(2, i))/Lambdas(1,i));
                    z(indexes) = den(i+1);
                elseif Lambdas(1,i) > 0 && Lambdas(1,i-1) > 0
                    indexes =  (x >= ((t- Lambdas(2, i))/Lambdas(1,i)) & x<= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                    z(indexes) = den(i);
                    indexes = (x <= (t- Lambdas(2, i))/Lambdas(1,i));
                    z(indexes) = den(i+1);
                elseif Lambdas(1,i) < 0 && Lambdas(1,i-1) > 0
                    indexes =  (x <= ((t- Lambdas(2, i))/Lambdas(1,i)) & x<= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                    z(indexes) = den(i);
                    indexes = (x >= (t- Lambdas(2, i))/Lambdas(1,i));
                    z(indexes) = den(i+1);
                elseif Lambdas(1,i) < 0 && Lambdas(1,i-1) < 0
                    indexes =  (x <= ((t- Lambdas(2, i-1))/Lambdas(1,i-1))) & (x <= (t- Lambdas(2, i))/Lambdas(1,i));
                    z(indexes) = den(i);
                    indexes = (x >= (t- Lambdas(2, i))/Lambdas(1,i));
                    z(indexes) = den(i+1);
                else
                    indexes =  (x >= ((t- Lambdas(2, i-1))/Lambdas(1,i-1))) & (x >= (t- Lambdas(2, i))/Lambdas(1,i));
                    z(indexes) = den(i);
                    indexes = (x <= (t- Lambdas(2, i))/Lambdas(1,i));
                    z(indexes) = den(i+1);

                end
            elseif i == 2
                if  Lambdas(1,i) > 0 && Lambdas(1,i-1) < 0
                    indexes =  (x >= ((t- Lambdas(2, i))/Lambdas(1,i)) & x>= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                    z(indexes) = den(i);
                    indexes = (x <= (t- Lambdas(2, i-1))/Lambdas(1,i-1));
                    z(indexes) = den(i-1);
                elseif Lambdas(1,i) > 0 && Lambdas(1,i-1) > 0
                    indexes =  (x >= ((t- Lambdas(2, i))/Lambdas(1,i)) & x<= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                    z(indexes) = den(i);
                    indexes = (x >= (t- Lambdas(2, i-1))/Lambdas(1,i-1));
                    z(indexes) = den(i-1);
                elseif Lambdas(1,i) < 0 && Lambdas(1,i-1) > 0
                    indexes =  (x <= ((t- Lambdas(2, i))/Lambdas(1,i)) & x<= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                    z(indexes) = den(i);
                    indexes = (x>= (t- Lambdas(2, i-1))/Lambdas(1,i-1));
                    z(indexes) = den(i-1);
                elseif Lambdas(1,i) == 0 || Lambdas(1,i-1) == 0
                    indexes =  (x <= ((t- Lambdas(2, i))/Lambdas(1,i)) & x>= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                    z(indexes) = den(i);
                    indexes = (x <= (t- Lambdas(2, i-1))/Lambdas(1,i-1));
                    z(indexes) = den(i-1);
                else
                    indexes =  (x <= ((t- Lambdas(2, i-1))/Lambdas(1,i-1))) & (x <= (t- Lambdas(2, i))/Lambdas(1,i));
                    z(indexes) = den(i);
                    indexes = (x <= (t- Lambdas(2, i-1))/Lambdas(1,i-1));
                    z(indexes) = den(i-1);
                end
            else
                if  Lambdas(1,i) > 0 && Lambdas(1,i-1) < 0
                    indexes =  (x >= ((t- Lambdas(2, i))/Lambdas(1,i)) & x>= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                elseif Lambdas(1,i) > 0 && Lambdas(1,i-1) > 0
                    indexes =  (x >= ((t- Lambdas(2, i))/Lambdas(1,i)) & x<= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                elseif Lambdas(1,i) < 0 && Lambdas(1,i-1) > 0
                    indexes =  (x <= ((t- Lambdas(2, i))/Lambdas(1,i)) & x<= ((t- Lambdas(2, i-1))/Lambdas(1,i-1)));
                else
                    indexes =  (x >= ((t- Lambdas(2, i-1))/Lambdas(1,i-1))) & (x <= (t- Lambdas(2, i))/Lambdas(1,i));
                end
                z(indexes) = den(i);
            end
            if width(Lambdas) == 2
                if Lambdas(1,1) > 0
                    indexes = (x >= (t- Lambdas(2, 1))/Lambdas(1,1));
                    z(indexes) = den(1);

                else
                    indexes = (x <= (t- Lambdas(2, 1))/Lambdas(1,1));
                    z(indexes) = den(1);

                end
            end

        end
        if width(Lambdas) == 1
            if Lambdas(1,1) > 0 %Plot density map
                indexes = (x >= (t- Lambdas(2, 1))/Lambdas(1,1));
                z(indexes) = den(1);
                indexes = (x <= (t- Lambdas(2, 1))/Lambdas(1,1));
                z(indexes) = den(2);
            else
                indexes = (x <= (t- Lambdas(2, 1))/Lambdas(1,1));
                z(indexes) = den(1);
                indexes = (x >= (t- Lambdas(2, 1))/Lambdas(1,1));
                z(indexes) = den(2);
            end
        end

        %Update Density List
        den = newden;

    end

    bigZ(int32(xIndexes),:) = z; %Fill in part of complete graph.

    end

if plotdecision == 1
    %Plot entire density map.
    figure("Color","White","Name","General LWR Solution")
    [bigX,bigT] = meshgrid(linspace(0,endDistance,endDistance/resolution),linspace(0,endTime,height(bigZ)));
    %subplot(1,2,1)
    surf(bigX,bigT,bigZ,  'EdgeColor', 'none')
    c = colorbar;
    c.Label.String = 'Density (veh/km)';
    c.Label.FontSize = 12;
    c.Label.FontWeight = "bold";
    c.Label.FontName = 'Times New Roman';
    view(2)
    xlabel('Distance, x (km)',"FontSize",12,"FontWeight","bold","FontName","Times New Roman")
    ylabel('Time t (minutes)',"FontSize",12,"FontWeight","bold","FontName","Times New Roman")
    xlim([0 endDistance])
    ylim([inputData(1,1) endTime])
    filename = sprintf('sim%s.jpg',fundamentalDiagram);
    exportgraphics(gcf,filename)

end

%% Find out predicted Sensor 2 densities.
outputDen = [];
for i = timeSteps(1):timeSteps(2) - timeSteps(1):endTime
    if i == 0
        extractIndex = 1;
    elseif i == endTime
        extractIndex = ((1/resolution)*i) - 1;
    else
        extractIndex = (1/resolution)*i;
    end
    density  = bigZ(extractIndex,width(bigZ));
    outputDen = [outputDen;[i-1,density]];

end
end

