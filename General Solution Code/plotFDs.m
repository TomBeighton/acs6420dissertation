function varargout = plotFDs(usage,method, sensor1,sensor2)
%plotFDs.m created by Tom Beighton for ACS6420 Advanced Project
%A function to plot Fundamental Diagrams and determine parameters from
%Epirical Fundamental Diagrams.

vf = 0.672; %km/min Traffic Free Flow Speed

switch method %Look at Fundamental Diagram Selected
    case "Greenshields"
        [rhoJam, ~,~,~ ]  = plotFDs("parameters","CustomOverall",sensor1,sensor2);
        %Call same function again to find parameters for tuning
        x = linspace(0,rhoJam,1000); 
        y = x.*vf.*(1-(x./rhoJam)); %Create vectors for cooridnates of FD plot
        
        %Plot all datapoints
        totalValues = [sensor1;sensor2];

        figure("Name",method,"Color","White")
        plot(totalValues(:,1),totalValues(:,2),'*') %Plot FD
        hold on
        plot(x,y);
        grid on

        xlabel("Density (Veh/Km)")
        ylabel("Flow (Veh/min)")

    case "Newell-Daganzo"

        [rhoJam, QMax,rhoCrit,~ ]  = plotFDs("parameters","CustomOverall",sensor1,sensor2);
        x = linspace(0,rhoJam,1000);
        y1 = (QMax/rhoCrit).*x(x<rhoCrit);
        y2 = (QMax/(rhoCrit - rhoJam)).*x(x>=rhoCrit) - (QMax*rhoJam)/(rhoCrit - rhoMax);

        y = [y1,y2];
        %Plot all datapoints
        totalValues = [sensor1;sensor2];
        figure("Name",method,"Color","White")
        plot(totalValues(:,1),totalValues(:,2),'*')
        hold on
        plot(x,y);
        grid on

        xlabel("Density (Veh/Km)")
        ylabel("Flow (Veh/min)")
    case "Greenberg"
        [rhoJam, ~,~,v0]  = plotFDs("parameters","CustomOverall",sensor1,sensor2);
        x = linspace(0,rhoJam,1000);
        y = x.*v0.*log(rhoJam./x);
        %Plot all datapoints
        totalValues = [sensor1;sensor2];
        figure("Name",method,"Color","White")
        plot(totalValues(:,1),totalValues(:,2),'*')
        hold on
        plot(x,y);
        grid on

        xlabel("Density (Veh/Km)")
        ylabel("Flow (Veh/min)")
        return

    case "Underwood"
        [rhoJam, ~,~,~ ]  = plotFDs("parameters","CustomOverall",sensor1,sensor2);
        x = linspace(0,rhoJam,1000);
        y =  x.*vf.*exp(-1.*(x/rhoJam));
        %Plot all datapoints
        totalValues = [sensor1;sensor2];
        figure("Name",method,"Color","White")
        plot(totalValues(:,1),totalValues(:,2),'*')
        hold on
        plot(x,y);
        grid on

        xlabel("Density (Veh/Km)")
        ylabel("Flow (Veh/min)")
        return

    case "California"
        [rhoJam, ~,~,v0 ]  = plotFDs("parameters","CustomOverall",sensor1,sensor2);
        x = linspace(0,rhoJam,1000);
        y = x.*v0.*((1./x)-(1/rhoJam));
        %Plot all datapoints
        totalValues = [sensor1;sensor2];
        figure("Name",method,"Color","White")
        plot(totalValues(:,1),totalValues(:,2),'*')
        hold on
        plot(x,y);
        grid on

        xlabel("Density (Veh/Km)")
        ylabel("Flow (Veh/min)")

        return
    case "CustomSensor1"
        %Plot both sensor 1
        totalValues = sensor1;
        figure("Name",method,"Color","White")
        plot(totalValues(:,1),totalValues(:,2),'*')
        hold on

        %Generate Trend line (Code created with MATLAB Curve Fitter App)
        % Set up fittype and options.
        ft = fittype( 'poly2' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Lower = [-Inf -Inf 0];
        opts.Upper = [Inf Inf 0];

        % Fit model to data.
        [fitresult,~] = fit( totalValues(:,1), totalValues(:,2), ft, opts );
        x = linspace(0,roots([fitresult.p1,fitresult.p2]),1000);
        y = fitresult.p1.*x.^2 + fitresult.p2.*x;

        if usage == "coeff"
            varargout{1} = fitresult.p1;
            varargout{2} = fitresult.p2;
        else
            % Plot fit with data.
            plot(x,y);
            grid on
            xlabel("Density (Veh/Km)")
            ylabel("Flow (Veh/min)")
        end


    case "CustomSensor2"
        %Plot  sensor 2 points
        totalValues = sensor2;
        figure("Name",method,"Color","White")
        plot(totalValues(:,1),totalValues(:,2),'*')
        hold on

        %Generate Trend line (Code created with MATLAB Curve Fitter App)
        % Set up fittype and options.
        ft = fittype( 'poly2' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Lower = [-Inf -Inf 0];
        opts.Upper = [Inf Inf 0];

        % Fit model to data.
        [fitresult,~] = fit( totalValues(:,1), totalValues(:,2), ft, opts );
        x = linspace(0,roots([fitresult.p1,fitresult.p2]),1000);
        y = fitresult.p1.*x.^2 + fitresult.p2.*x;

        if usage == "coeff"
            varargout{1} = fitresult.p1;
            varargout{2} = fitresult.p2;
        else
            % Plot fit with data.
            plot(x,y);
            grid on
            xlabel("Density (Veh/Km)")
            ylabel("Flow (Veh/min)")
        end


    case "CustomOverall"
        %Collate both sensor 1 and sensor 2 points
        totalValues = [sensor1;sensor2];


        %Generate Trend line (Code created with MATLAB Curve Fitter App)
        % Set up fittype and options.
        ft = fittype( 'poly2' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Lower = [-Inf -Inf 0];
        opts.Upper = [Inf Inf 0];

        % Fit model to data.
        [fitresult,~] = fit( totalValues(:,1), totalValues(:,2), ft, opts );
        x = linspace(0,roots([fitresult.p1,fitresult.p2]),1000);
        y = fitresult.p1.*x.^2 + fitresult.p2.*x;


        if usage == "parameters"
            %Figure out QMax, RhoCrit, RhoJam
            rhoJam = roots([fitresult.p1,fitresult.p2]);
            QMax = max(y);
            rhoCrit = roots([fitresult.p1,fitresult.p2,-QMax]);
            rhoCrit = rhoCrit(1);
            v0 = QMax / rhoCrit;
            %Return Parameters
            varargout{1} = rhoJam;
            varargout{2} = QMax;
            varargout{3} = rhoCrit;
            varargout{4} = v0;
        elseif usage == "plot"
            % Plot fit with data.
            figure("Name",method,"Color","White")
            plot(totalValues(:,1),totalValues(:,2),'*')
            hold on
            plot(x,y);
            grid on

            xlabel("Density (Veh/Km)")
            ylabel("Flow (Veh/min)")
        elseif usage == "coeff"
            varargout{1} = fitresult.p1;
            varargout{2} = fitresult.p2;
        end

end

end
