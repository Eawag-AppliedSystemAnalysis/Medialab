function plotFred(FredContribution, speciesName, plBC, simValues, time, isSolidSpecies, units, porosity, rho, o2Input)

% plotFred plots the flux of reduced substances out of the sediment, Fred, 
% as a function of time, along with OM and O2 BCs.
%
%
% input arguments are:
%   FredContribution    - array specifying contribution to Fred for each species
%   speciesName         - cell array containing strings of species
%                         names (or name of reactions etc.). Name of species end with (s) or (aq)
%                         to distinguish between solid or
%                         aqeouse species, respectively.
%   plBC                - function handles for fluxes at the upper boundary
%   simValues           - stores the simulated concentrations (or other variable) of all the 
%                         species at all depths and times. Its size is
%                         {numSpecies*1} cells where each cell stores
%                         concentrations of the species at nX depths and nT 
%                         times i.e. [nT*nX] size matrices.
%   time                - temporal interval in [yr] with nT discretized
%                         steps
%   isSolidSpecies      - specifies which species are solids and solutes
%   units               - string or cell array of strings representing the units that should be used
%                         to label the x-axis. If two are given, the first
%                         is assumed to be for solid and the second for
%                         aqueous species
%    porosity           - porosity of the sediment
%    rho                - dry density of the sediment         
%    o2Input            - function handle taking time as an argument and
%                         returning the upper BC of O2 as output


% Initialization
nSpecies = length(speciesName);   
nTimes = min(length(time), size(simValues{1},1)); % To work even if the solver aborted before finishing
time = time(1:nTimes);
clrs = jet(length(speciesName));
lineStyles = {':','-','--','-.'};
figure

% Extract concentrations at upper boundary
for iSpecies=1:nSpecies
    ul(iSpecies,:) = simValues{iSpecies}(:,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***************************CALCULATE FLUXES******************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSpecies=1:nSpecies
    % fluxes at upper and lower boundary for all species
    for iTime=1:nTimes
        if isSolidSpecies(iSpecies)
            fluxl{iSpecies}(iTime) = plBC{iSpecies}(ul(:,iTime),time(iTime))*rho*(1-porosity);
        else
            fluxl{iSpecies}(iTime) = plBC{iSpecies}(ul(:,iTime),time(iTime))*porosity*(1.02*porosity^-1.81);
        end
    end
   
    
    % fluxes out of sediment
    soluteFlux = -FredContribution(iSpecies)*(fluxl{iSpecies})*(~isSolidSpecies(iSpecies));
    solidFlux = zeros(size(soluteFlux));
    fluxOut(iSpecies,:) = max(0,solidFlux+soluteFlux);
end

fluxOut = fluxOut/10^6*10000/365.25*16*2; % To change units from mumol/cm2/yr to g O2/m2/d


om1Flux = fluxl{12};
o2Flux = fluxl{1}/10^6*10000/365.25*16*2;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%****************************PLOT FRED************************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,1,1)
set(gcf,'Position',get(0,'ScreenSize'));                          %maximize the figure
set(gcf,'PaperPosition',[3,5,24,11]);                             %set the dimensions on paper when the plots are printed
hold on;
plot(time,sum(fluxOut),'r')
legendEntries = {'Total F_{red}'};
for iSpecies = 1:nSpecies
    if FredContribution(iSpecies)>0
        plot(time,fluxOut(iSpecies,:),'color',clrs(iSpecies,:),'LineStyle',lineStyles{mod(iSpecies,4)+1})
        legendEntries = {legendEntries{:}, speciesName{iSpecies}};
    end
end
plot(time,o2Flux,'b')
plot(time,sum(fluxOut)+o2Flux,'k')
legendEntries = {legendEntries{:}, 'O2 consumption','O2 + F_{red}'};

legend(legendEntries,'Location','Best')
title('F_{red}', 'FontSize',16,'FontWeight','bold');
xlabel('time [yr]','FontSize',14);
ylabel(units,'FontSize',14);
% ylim([0,1.05*max(max(sum(fluxIn)),max(sum(fluxOut)))])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************PLOT O2 BC********************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
subplot(3,1,2)
hold on;
for iTime = 1:length(time)
    o2BC(iTime) = o2Input(time(iTime));
end
plot(time,o2BC)
title('Boundary Condition O2', 'FontSize',16,'FontWeight','bold');
xlabel('time [yr]','FontSize',14);
ylabel('[\mumol/cm^3]','FontSize',14);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******************************PLOT OM1 BC*******************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
subplot(3,1,3)
hold on;
plot(time,om1Flux)
title('Boundary Condition OM1', 'FontSize',16,'FontWeight','bold');
xlabel('time [yr]','FontSize',14);
ylabel('[\mumol/cm^2/yr]','FontSize',14);

SumFred=sum(fluxOut)
SumO2=o2Flux
SumAhm=sum(fluxOut)+o2Flux
save('AHM_Fred_O2.mat','SumFred','SumO2','SumAhm')
