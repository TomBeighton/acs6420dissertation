%ACS6420 Advanced Project: Tom Beighton 180160767
%FDs.m : Function allowing the generation of a Lambda variable based on
%        density inputs from GeneralSolution.m and a chosen method.


function [lambda] = Fds(rhoMinus,rhoPlus,method, vf, v0, rhoMax,QMax,rhoCrit,p1,p2)
    
    switch method
        case "TrafficFreeFlow"
            QPlus = vf*rhoPlus;
            QMinus = vf*rhoMinus;
            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);
        case "Greenshields"
            QPlus = rhoPlus*vf*(1 - (rhoPlus/rhoMax));
            QMinus = rhoMinus*vf*(1 - (rhoMinus/rhoMax));
            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);


            return 
        case "Newell-Daganzo"
            if rhoPlus > rhoCrit
                QPlus = (QMax/(rhoCrit - rhoMax))*rhoPlus - (QMax*rhoMax)/(rhoCrit - rhoMax);
            else
                QPlus = (QMax /rhoCrit);
            end

            if rhoMinus > rhoCrit
                QMinus = (QMax/(rhoCrit - rhoMax))*rhoMinus - (QMax*rhoMax)/(rhoCrit - rhoMax);
            else
                QMinus = (QMax /rhoCrit);
            end
            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);

            return
    
        case "Greenberg"
            QPlus = rhoPlus*v0*log(rhoMax/rhoPlus);
            QMinus = rhoMinus*v0*log(rhoMax/rhoMinus);
            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);
            return

        case "Underwood"
            QPlus = rhoPlus*vf*exp(-1*(rhoPlus/rhoMax));
            QMinus = rhoMinus*vf*exp(-1*(rhoPlus/rhoMax));
            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);
            return

        case "California"

            QPlus = rhoPlus*v0*((1/rhoPlus)-(1/rhoMax));
            QMinus = rhoMinus*v0*((1/rhoMinus)-(1/rhoMax));
            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);
            return
        case "CustomSensor1"
            %Values drawn from Fundamental Diagram Experiments
            QPlus = p1*(rhoPlus)^2 + p2*rhoPlus;
            QMinus =p1*(rhoMinus)^2 +p2*rhoMinus;

            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);
            return
        case "CustomSensor2"
            QPlus = p1*(rhoPlus)^2 + p2*rhoPlus;
            QMinus =p1*(rhoMinus)^2 +p2*rhoMinus;

            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);
            return
        case "CustomOverall"
            QPlus = p1*(rhoPlus)^2 + p2*rhoPlus;
            QMinus =p1*(rhoMinus)^2 +p2*rhoMinus;

            lambda = (QPlus - QMinus)/(rhoPlus - rhoMinus);
            return
   end

end


