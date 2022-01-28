function plotWithMeasurements(speciesName, isSolidSpecies, simValues, time, depth, fieldDataFile, dataName, speciesToPlot)

% plotWithMeasurements plots the measured and modeled concentration 
% profiles of species for several points in time
%
%
% input arguments are:
%   speciesName         - cell array containing strings of species
%                         names. Name of species end with (s) or (aq)
%                         to distinguish between solid or 
%                         aqeouse species, respectively. Size is nSpecies
%   isSolidSpecies      - 1 if species is solid, 0 if it is solute (is
%                         relevant for correct labels). Size is nSpecies
%   simValues           - stores the simulated concentrations of all the 
%                         species at all depths and times. Its size is
%                         {nSpecies*1} cells where each cell stores
%                         concentrations of the species at nX depths and nT 
%                         times i.e. [nT*nX] size metrices.
%   time                - time in [yr] of simValues
%   depth               - spatial interval in [cm] with nX descritized 
%                         depths
%   fieldDataFile       - MATLAB file with measured data included
%                         in a specific format.
%   dataName            - string variable of species name to be plotted  
%   speciesToPlot       - components of dataName.
% ------------------------------------------------------------------------%


% Initialization
load(fieldDataFile, 'profiles');
nSpeciesToPlot = length(speciesToPlot);
clrs = jet(length(profiles));

for s = 1:nSpeciesToPlot
    speciesToPlotNoPhase{s} = regexprep(speciesToPlot{s},'\(.+\)',''); % names of species without the appended phase information
end



% Generate the plot(s)
figure; 
set(gcf,'Position',get(0,'ScreenSize'));                                  % maximize the figure
set(gcf,'PaperPosition',[3,5,24,11]);                                     % set the dimensions on paper when the plots are printed
set(gca,'FontSize',18);                                                    % set the font size

for iSpecies = 1:nSpeciesToPlot
    subplot(2,2,iSpecies);
    hold on
    
    legendEntries = {};
    index = find(strcmp(speciesName,speciesToPlot{iSpecies}) ~= 0);
    for iProfile = 1:length(profiles)
        [~, step] = min(abs(time-(floor(max(time))-1+profiles(iProfile).time)));
        y = simValues{index}(step,:);
        plot(y,depth,'LineWidth',1.5,'color',clrs(iProfile,:));
        scatter(profiles(iProfile).(speciesToPlotNoPhase{iSpecies}).concentration, profiles(iProfile).(speciesToPlotNoPhase{iSpecies}).depth, [],clrs(iProfile,:)); 
        legendEntries = {legendEntries{:},sprintf('Model (%.2f yr)',profiles(iProfile).time),'Data'};
    end
    
    set(gca,'YDir','reverse');                                                % reverse the y-axis, since depth should increase from top to bottom
    ylim([0 depth(end)])
    
    % Label the plot
    if isSolidSpecies
        xlabel('\mumol/g','FontSize',18);                                       % x axis label for solid
    else
        xlabel('\mumol/cm^3','FontSize',18);                                    % x axis label for solute
    end
    
    ylabel('Depth (cm)','FontSize',18);                                        % label of y-axis
    
    title(char(speciesToPlot(iSpecies)), 'FontSize',18,'FontWeight','bold');                 % title of the plots is the name of the species

    legend(legendEntries,'Location','Best') 
end