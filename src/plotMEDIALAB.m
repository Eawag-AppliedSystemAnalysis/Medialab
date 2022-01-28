function plotMEDIALAB(speciesName, simValues,depth, fieldDataFile, dataName, speciesToPlot)

% plotMEDIALAB plots the measured and modeled concentration profiles of 
% species
%
%
% input arguments are:
%   speciesName         - cell array containing strings of species
%                         names. Name of species end with (s) or (aq)
%                         to distinguish between solid or
%                         aqeouse species, respectively.
%   simValues           - stores the simulated concentrations of all the 
%                         species at all depths and times. Its size is
%                         {numSpecies*1} cells where each cell stores
%                         concentrations of the species at nX depths and nT 
%                         times i.e. [nT*nX] size metrices.
%   depth               - spatial interval in [cm] with nX discretized 
%                         depths
%   fieldDataFile       - excel file with measured data included
%                         in a specific format.
%   dataName            - string variable of species name to be plotted  
%   speciesToPlot       - components of dataName.


x            = depth;
[m,n]        = size(simValues{1,1}(:,:));
simValuesNew = zeros(m,n);

if (isempty(strfind(dataName,'(aq)'))); 
    speciesPhase = 1;
else
    speciesPhase = 0;
end
   
if (length(speciesToPlot) > 1)
    for i = 1:length(speciesToPlot)
        titleName = dataName;
        index = find(strcmp(speciesName,speciesToPlot{i}) ~= 0);
        if ((strcmp(speciesToPlot{i},'om1(s)')~= 0)&&(strcmp(dataName,'TP(s)')~= 0))
            simValuesNew = simValuesNew + simValues{index,1}(:,:)/106;     % 106 is the C/P ratio of OM1
        elseif ((strcmp(speciesToPlot{i},'om2(s)')~= 0)&&(strcmp(dataName,'TP(s)')~= 0))
            simValuesNew = simValuesNew + simValues{index,1}(:,:)/106;     % 106 is C/P ratio of OM2
        elseif ((strcmp(speciesToPlot{i},'om3(s)')~= 0)&&(strcmp(dataName,'TP(s)')~= 0))
            simValuesNew = simValuesNew + simValues{index,1}(:,:)/356;     % 356 is C/P ratio of OM3
        elseif ((strcmp(speciesToPlot{i},'vivianite(s)')~= 0)&&(strcmp(dataName,'TP(s)')~= 0))
            simValuesNew = simValuesNew + simValues{index,1}(:,:)*2;       % 2 is numner of P moles in vivianite
        elseif ((strcmp(speciesToPlot{i},'vivianite(s)')~= 0)&&(strcmp(dataName,'TFe(s)')~= 0))
            simValuesNew = simValuesNew + simValues{index,1}(:,:)*3;       % 3 is number of Fe moles in vivianite
        else
            simValuesNew = simValuesNew + simValues{index,1}(:,:);
        end
    end
else
        titleName = speciesToPlot;
        index = find(strcmp(speciesName,speciesToPlot) ~= 0);
        if ((strcmp(speciesToPlot,'om1(s)')~= 0)&&(strcmp(dataName,'TP(s)')~= 0))
            simValuesNew = simValues{index,1}(:,:)/106;
        elseif ((strcmp(speciesToPlot,'om2(s)')~= 0)&&(strcmp(dataName,'TP(s)')~= 0))
            simValuesNew = simValues{index,1}(:,:)/106;
        elseif ((strcmp(speciesToPlot,'om3(s)')~= 0)&&(strcmp(dataName,'TP(s)')~= 0))
            simValuesNew = simValues{index,1}(:,:)/356;
        elseif ((strcmp(speciesToPlot,'vivianite(s)')~= 0)&&(strcmp(dataName,'TP(s)')~= 0))
            simValuesNew = simValues{index,1}(:,:)*2;
        elseif ((strcmp(speciesToPlot,'vivianite(s)')~= 0)&&(strcmp(dataName,'TFe(s)')~= 0))
            simValuesNew = simValues{index,1}(:,:)*3;    
        else
            simValuesNew = simValues{index,1}(:,:);
        end

end


 y = simValuesNew(1,:);
 plot(y,x,'LineWidth',1.5,'color','b');
 hold on;
 
 y = simValuesNew(floor(end/2),:);
 plot(y,x,'LineWidth',1.5,'color','m');
 
 y = simValuesNew(end,:);
 plot(y,x,'LineWidth',1.5,'color','r');
 
 set(gca,'YDir','reverse');                                                % reverse the y-axis, since depth should increase from top to bottom
 set(gcf,'Position',get(0,'ScreenSize'));                                  % maximize the figure
 set(gcf,'PaperPosition',[3,5,24,11]);                                     % set the dimensions on paper when the plots are printed
 set(gca,'FontSize',18);                                                    % set the font size

    data = xlsread(fieldDataFile, char(dataName));                         % stores the data in excel sheet in the cell matrix Data
    if (size(data,2) >= 2)
        scatter(data(:,2), data(:,1) );                                    % plots a scatter diagram of conc. vs. depth                  
    end

% xl = xlim;
% if xl(1)<0
%    xlim([0,xl(2)]);                                                        % set the x-axis to start from 0 (if it's negative)
% end

% xlim([0 7000])
ylim([0 depth(end)])
if speciesPhase==1
   xlabel('\mumol/g','FontSize',18);                                       % x axis label for solid
 else
   xlabel('\mumol/cm^3','FontSize',18);                                    % x axis label for solute
end

ylabel('Depth (cm)','FontSize',18);                                        % label of y-axis

legend('IC','1/2 time','end (eq)','meas','Location','Best')

 title(char(titleName), 'FontSize',18,'FontWeight','bold');                 % title of the plots is the name of the species
