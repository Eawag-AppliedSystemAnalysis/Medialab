function plotAgainstDepth(titleText, speciesName, simValues, depth, units, useLogscale)

% plotAgainstDepth plots modeled profiles of species at all times specified in
% userMEDIALAB
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
%   depth               - spatial interval in [cm] with nX descritized 
%                         depths
%   units               - string or cell array of strings representing the units that should be used
%                         to label the x-axis. If two are given, the first
%                         is assumed to be for solid and the second for
%                         aqueous species
%   useLogscale            - boolean determining whether normal values are
%                         plotted or the log10
%                         


% Initialization
N = length(simValues);
clrs = jet(size(simValues{1},1));
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
        for j = 1:size(simValues{i},1)
            x = simValues{i}(j,:);
            if useLogscale
                x = log(x)/log(10).*(x>0);
            end
            plot(x,depth,'color',clrs(j,:));
        end
        set(gca,'YDir','reverse');                                        %reverse the y-axis, since depth should increase from top to bottom
        
        
        % Define the x- and y-limits
        xmax=max(max(simValues{i}));
        if(xmax==0);xmax=0.5;end
        xmin = 0;
        if useLogscale
            xmax = log(xmax)/log(10);
            minTemp = min(min(simValues{i}));
            if minTemp>0
                xmin = log(minTemp)/log(10);
            else
                if minTemp==0
                    xmin = xmax-3;
                else
                    xmax = 0.5;
                    xmin = min(xmax-3,log(max(max(-simValues{i})))/log(10));
                end
            end
        end
        try
        xlim([xmin xmax]);
        catch
            
        end
        
        
        % Label the figure
        title(char(speciesName{i}), 'FontSize',16,'FontWeight','bold');
        if mod(i,16)==1
            legend('IC','Location','Best')
        end
        if mod(i,4)==1
            ylabel('depth [cm]','FontSize',14);
        end
        if length(units)==1
            if (i>=(N-4))||(mod(i,16)>=13)
                xlabel(units,'FontSize',14);
            end
        else
            if (isempty(strfind(speciesName{i},'(aq)')));
                xlabel(units{1},'FontSize',14);
            else
                xlabel(units{2},'FontSize',14);
            end
        end
    end
    
    % Plot the title in versions where this works
    if strcmp(verMatlab.Version,'7.10')
        ha = axes('Parent',hp,'Visible','off');
        ht = text('Parent',ha, 'String',titleText,'FontSize',18,'FontWeight','bold','Position',[0.5 1.05]);
    end
end











