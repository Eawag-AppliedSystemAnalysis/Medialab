function plotAgainstTime(titleText, speciesName, simValues, time, depth, iDepthToPlot, units, useLogscale)

% plotAgainstTime plots how concentrations, rates, etc. change over time at
% specified depths
%
%
% input arguments are:
%   titleText           - title of the figure
%   speciesName         - cell array containing strings of species
%                         names (or name of reactions etc.). Name of species end with (s) or (aq)
%                         to distinguish between solid or
%                         aqeouse species, respectively.
%   simValues           - stores the simulated concentrations (or other variable) of all the 
%                         species at all depths and times. Its size is
%                         {numSpecies*1} cells where each cell stores
%                         concentrations of the species at nX depths and nT 
%                         times i.e. [nT*nX] size matrices.
%   time                - temporal interval in [yr] with nT discretized
%                         steps
%   depth               - spatial interval in [cm] with nX descritized 
%                         depths
%   iDepthToPlot        - array of indices, specifying which depths to plot
%   units               - string or cell array of strings representing the units that should be used
%                         to label the x-axis. If two are given, the first
%                         is assumed to be for solid and the second for
%                         aqueous species
%   useLogscale            - boolean determining whether normal values are
%                         plotted or the log10
%                         


% Initialization
clrToPlot = {'b','r','k','c','m','g','y'};
depthLabel = {};
for j = 1:length(iDepthToPlot)
    depthLabel{j} = sprintf('depth %0.2f cm',depth(iDepthToPlot(j)));
end
N = length(simValues);
if useLogscale
    for i=1:N
        speciesName{i} = sprintf('log(%s)',speciesName{i});
    end
    for i = 1:length(units)
        units{i} = sprintf('log_{10}(%s)',units{i});
    end
end


% Generate as many figures as necessary with up to 16 panels
verMatlab = ver('MATLAB'); % check the version of Matlab
for n=1:ceil(N/16)
    figure;
    % Prepare uipanel to plot the title of the figure in the end
    if strcmp(verMatlab.Version,'7.10')
        hp = uipanel;   
        set(hp,'BackgroundColor','none') %make uipanel invisible
    end
    set(gcf,'Position',get(0,'ScreenSize'));                          %maximize the figure
    set(gcf,'PaperPosition',[3,5,24,11]);                             %set the dimensions on paper when the plots are printed
    for i=(1+(n-1)*16):min(n*16,N)
        subplot(4,4,i-(n-1)*16)
        set(gca,'FontSize',12);                                            %set the font size
        hold on;
        for j = 1:length(iDepthToPlot)
            y = simValues{i}(:,iDepthToPlot(j));
            if useLogscale
                y = log(y)/log(10).*(y>0);
            end
            plot(time,y,'color',clrToPlot{j});
        end
        
        
        % Define the x- and y-limits
        ymin = min(min(simValues{i}(:,iDepthToPlot)));
        ymax = max(max(simValues{i}(:,iDepthToPlot)));
        if useLogscale
            ymin = log(ymin)/log(10);
            ymax = log(ymax)/log(10);
            if isinf(ymin)
                ymin = ymax-3;
            end
        end
        ylim([min(0, ymin) max(0.00000001, ymax)])
        
        
        % Label the figure
        title(char(speciesName{i}), 'FontSize',16,'FontWeight','bold');
        if mod(i,16)==1
            legend(depthLabel,'Location','Best')
        end
        if (i>=(N-4))||(mod(i,16)>=13)
            xlabel('time [yr]','FontSize',14);
        end
        if length(units)==1
            if mod(i,4)==1
                ylabel(units,'FontSize',14);
            end
        else
            if (isempty(strfind(speciesName{i},'(aq)')));
                ylabel(units{1},'FontSize',14);
            else
                ylabel(units{2},'FontSize',14);
            end
        end
    end
    
    % Plot the title in versions where this works
    if strcmp(verMatlab.Version,'7.10')
        ha = axes('Parent',hp,'Visible','off');
        ht = text('Parent',ha, 'String',titleText,'FontSize',18,'FontWeight','bold','Position',[0.5 1.05]);
    end
end











