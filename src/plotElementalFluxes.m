function plotElementalFluxes(elementName, elementContent, speciesName, plBC, prBC, simValues, time, depth, isSolidSpecies, units, porosity, rho)

% plotElementalFluxes plots the elemental fluxes in and out of the modeled
% sediment slab and calculates the mass balance of these fluxes.
% Additionally, the temporal change of the inventory of each element in 
% the modeled slab is calculated, which allows to check whether the mass
% balance is correct.
%
%
% input arguments are:
%   elementName         - Name of the element e.g. 'C'
%   elementContent      - array specifying how many atoms of the element
%                         are present in each species
%   speciesName         - cell array containing strings of species
%                         names (or name of reactions etc.). Name of species end with (s) or (aq)
%                         to distinguish between solid or
%                         aqeouse species, respectively.
%   plBC, prBC          - function handles for fluxes at boundaries
%   simValues           - stores the simulated concentrations (or other variable) of all the 
%                         species at all depths and times. Its size is
%                         {numSpecies*1} cells where each cell stores
%                         concentrations of the species at nX depths and nT 
%                         times i.e. [nT*nX] size matrices.
%   time                - temporal interval in [yr] with nT discretized
%                         steps
%   depth               - spatial interval in [cm] with nX descritized 
%                         depths
%   isSolidSpecies      - specifies which species are solids and solutes
%   units               - string or cell array of strings representing the units that should be used
%                         to label the x-axis. If two are given, the first
%                         is assumed to be for solid and the second for
%                         aqueous species
%    porosity           - porosity of the sediment
%    rho                - dry density of the sediment                      


% Initialization
nSpecies = length(speciesName);   
nTimes = min(length(time), size(simValues{1},1)); % To work even if the solver aborted before finishing
time = time(1:nTimes);
clrs = jet(length(speciesName));
lineStyles = {':','-','--','-.'};

% Extract concentrations at upper/lower boundary
for iSpecies=1:nSpecies
    ul(iSpecies,:) = simValues{iSpecies}(:,1);
    ur(iSpecies,:) = simValues{iSpecies}(:,length(depth));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***************************CALCULATE FLUXES******************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSpecies=1:nSpecies
    % fluxes at upper and lower boundary for all species
    for iTime=1:nTimes
        if isSolidSpecies(iSpecies)
            fluxl{iSpecies}(iTime) = plBC{iSpecies}(ul(:,iTime),time(iTime))*rho*(1-porosity);
            fluxr{iSpecies}(iTime) = prBC{iSpecies}(ur(:,iTime),time(iTime))*rho*(1-porosity);
        else
            fluxl{iSpecies}(iTime) = plBC{iSpecies}(ul(:,iTime),time(iTime))*porosity;
            fluxr{iSpecies}(iTime) = prBC{iSpecies}(ur(:,iTime),time(iTime))*porosity;
        end
    end
    
    % fluxes into sediment
    fluxIn(iSpecies,:) = max(0,elementContent(iSpecies)*fluxl{iSpecies});
    
    % fluxes out of sediment
    soluteFluxTop = -elementContent(iSpecies)*fluxl{iSpecies}*(~isSolidSpecies(iSpecies));
    soluteFluxBottom = elementContent(iSpecies)*fluxr{iSpecies}*(~isSolidSpecies(iSpecies));
    solidFlux = elementContent(iSpecies)*fluxr{iSpecies}*isSolidSpecies(iSpecies);
    fluxOutTop(iSpecies,:) = max(0,soluteFluxTop);
    fluxOutBottom(iSpecies,:) = max(0,solidFlux+soluteFluxBottom);
end

% Flux balance
totalFlux = sum(fluxIn)-sum(fluxOutTop)-sum(fluxOutBottom);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***********************CALCULATE MASS IN SEDIMENT************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSpecies=1:nSpecies
    for iTime=1:nTimes
        if isSolidSpecies(iSpecies)
            elementMass(iSpecies,iTime) = elementContent(iSpecies)*trapz(depth,simValues{iSpecies}(iTime,:))*rho*(1-porosity);
        else
            elementMass(iSpecies,iTime) = elementContent(iSpecies)*trapz(depth,simValues{iSpecies}(iTime,:))*porosity;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***********************FLUXES INTO SEDIMENT******************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;     
set(gcf,'Position',get(0,'ScreenSize'));                          %maximize the figure
set(gcf,'PaperPosition',[3,5,24,11]);                             %set the dimensions on paper when the plots are printed
    
subplot(2,3,1)
hold on;
plot(time,sum(fluxIn),'k')
legendEntries = {'Total'};
for iSpecies = 1:nSpecies
    if elementContent(iSpecies)>0
        plot(time,fluxIn(iSpecies,:),'color',clrs(iSpecies,:),'LineStyle',lineStyles{mod(iSpecies,4)+1})
        legendEntries = {legendEntries{:}, speciesName{iSpecies}};
    end
end
legend(legendEntries,'Location','Best')
title(sprintf('%s fluxes into sediment',elementName), 'FontSize',16,'FontWeight','bold');
xlabel('time [yr]','FontSize',14);
ylabel(units,'FontSize',14);
ylim([0,1.05*max(max(sum(fluxIn)),max(sum(fluxOutTop+fluxOutBottom)))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**********************FLUXES OUT OF SEDIMENT*****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2)
hold on;
plot(time,sum(fluxOutTop),'k')
legendEntries = {'Total'};
for iSpecies = 1:nSpecies
    if (elementContent(iSpecies)>0)&&(~isSolidSpecies(iSpecies))
        plot(time,fluxOutTop(iSpecies,:),'color',clrs(iSpecies,:),'LineStyle',lineStyles{mod(iSpecies,4)+1})
        legendEntries = {legendEntries{:}, speciesName{iSpecies}};
    end
end
legend(legendEntries,'Location','Best')
title(sprintf('%s fluxes out of sediment at top',elementName), 'FontSize',16,'FontWeight','bold');
xlabel('time [yr]','FontSize',14);
ylabel(units,'FontSize',14);
ylim([0,1.05*max(max(sum(fluxIn)),max(sum(fluxOutTop+fluxOutBottom)))])


subplot(2,3,4)
hold on;
plot(time,sum(fluxOutBottom),'k')
legendEntries = {'Total'};
for iSpecies = 1:nSpecies
    if (elementContent(iSpecies)>0)
        plot(time,fluxOutBottom(iSpecies,:),'color',clrs(iSpecies,:),'LineStyle',lineStyles{mod(iSpecies,4)+1})
        legendEntries = {legendEntries{:}, speciesName{iSpecies}};
    end
end
legend(legendEntries,'Location','Best')
title(sprintf('%s fluxes out of sediment at bottom',elementName), 'FontSize',16,'FontWeight','bold');
xlabel('time [yr]','FontSize',14);
ylabel(units,'FontSize',14);
ylim([0,1.05*max(max(sum(fluxIn)),max(sum(fluxOutTop+fluxOutBottom)))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************** SOLID MASS IN SEDIMENT****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
subplot(2,3,5)
hold on;
% plot(time,isSolidSpecies*elementMass,'k')
plot(time,ones(size(isSolidSpecies))*elementMass,'k')
legendEntries = {'Total solids'};
for iSpecies = 1:nSpecies
%     if (isSolidSpecies(iSpecies)&&elementContent(iSpecies)>0)
    if (elementContent(iSpecies)>0)
        plot(time,elementMass(iSpecies,:),'color',clrs(iSpecies,:),'LineStyle',lineStyles{mod(iSpecies,4)+1})
        legendEntries = {legendEntries{:}, speciesName{iSpecies}};
    end
end
legend(legendEntries,'Location','Best')
title(sprintf('%s mass in sediment',elementName), 'FontSize',16,'FontWeight','bold');
xlabel('time [yr]','FontSize',14);
ylabel(units(1:end-3),'FontSize',14);
ylim([0,1.05*max(ones(size(isSolidSpecies))*elementMass)])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**********************SOLUTE MASS IN SEDIMENT****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% set(gcf,'Position',get(0,'ScreenSize'));                          %maximize the figure
% set(gcf,'PaperPosition',[3,5,24,11]);                             %set the dimensions on paper when the plots are printed
%     
% subplot(2,3,5)
% hold on;
% plot(time,~isSolidSpecies*elementMass,'k')
% legendEntries = {'Total solutes'};
% for iSpecies = 1:nSpecies
%     if (~isSolidSpecies(iSpecies)&&elementContent(iSpecies)>0)
%         plot(time,elementMass(iSpecies,:),'color',clrs(iSpecies,:),'LineStyle',lineStyles{mod(iSpecies,4)+1})
%         legendEntries = {legendEntries{:}, speciesName{iSpecies}};
%     end
% end
% legend(legendEntries,'Location','Best')
% title(sprintf('%s mass in sediment (solutes)',elementName), 'FontSize',16,'FontWeight','bold');
% xlabel('time [yr]','FontSize',14);
% ylabel(units(1:end-3),'FontSize',14);
% ylim([0,1.05*max(max(sum(elementMass)),max(sum(elementMass)))])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***********************TOTAL MASS BALANCE********************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,3,3)
hold on;
legendEntries = {'Flux Balance','Sediment Mass Balance','Total Mass Balance'};

plot(time,totalFlux,'b')
plot(time,gradient(sum(elementMass),time),'r')
plot(time,totalFlux-gradient(sum(elementMass),time),'k')
legend(legendEntries,'Location','Best')
title(sprintf('%s Balance',elementName), 'FontSize',16,'FontWeight','bold');
xlabel('time [yr]','FontSize',14);
ylabel(units,'FontSize',14);

