function plotGradients(simValues, depth, time, selTime, speciesName)

% plotGradients plots the relative concentration gradient (Delta C/C) for 
% each species or reaction.
%
%
% input arguments are:
%   simValues   	- The model output 
%   depth           - Vector of depths at which simValues contains outputs
%   time            - vector containing all times
%   selTime         - Vector of indices selecting for which times to
%                     generate the plot. If it is empty, the code iterates
%                     through all times to identify the time at which the
%                     gradient was highest (for each depth) and plots this.
%                     plotted. Dimensions [length(X)*nComponents*nPlots].
%   speciesName     - names of the species
% ------------------------------------------------------------------------%



clrs = jet(length(speciesName));

% If selTime is empty, iterate through all time to identify the highest
% relative gradient at each depth
if isempty(selTime)
    
    figure;
    set(gcf,'Position',get(0,'ScreenSize'));                                  % maximize the figure
    set(gcf,'PaperPosition',[3,5,24,11]);                                     % set the dimensions on paper when the plots are printed
    hold on
    
    for iSpecies = 1:length(speciesName)
        for iTimes = 1:size(simValues{1},1)
            absGrad = abs(diff(simValues{iSpecies}(iTimes,:)));
            relGrad(iTimes,:) = absGrad.*2./(simValues{iSpecies}(iTimes,1:end-1)+simValues{iSpecies}(iTimes,2:end));
        end
        if mod(iSpecies,4)==1
            plot(max(relGrad),depth(1:end-1),'color',clrs(iSpecies,:))
        else if mod(iSpecies,4)==2
                plot(max(relGrad),depth(1:end-1),'color',clrs(iSpecies,:),'LineStyle','--')
            else if mod(iSpecies,4)==3
                    plot(max(relGrad),depth(1:end-1),'color',clrs(iSpecies,:),'LineStyle','-.')
                else
                    plot(max(relGrad),depth(1:end-1),'color',clrs(iSpecies,:),'LineStyle',':')
                end
            end
        end
    end
    
    set(gca,'YDir','reverse');                                                % reverse the y-axis, since depth should increase from top to bottom
    set(gca,'FontSize',18);                                                    % set the font size
    legend(speciesName,'NumColumns',2,'Location','Best','FontSize',12)
    xlabel('\Delta C/C','FontSize',16)
    ylabel('Depth [cm]','FontSize',16)
    title('Relative Concentration Gradient (Max)','FontSize',16,'FontWeight','bold')
    
    
% plot relative gradient at selected times
else
    for iTimes = selTime
        
        figure;
        hold on
        
        for iSpecies = 1:length(speciesName)
            absGrad = abs(diff(simValues{iSpecies}(iTimes,:)));
            relGrad = absGrad.*2./(simValues{iSpecies}(iTimes,1:end-1)+simValues{iSpecies}(iTimes,2:end));            
            
            if mod(iSpecies,4)==1
                plot(relGrad,depth(1:end-1),'color',clrs(iSpecies,:))
            else if mod(iSpecies,4)==2
                    plot(relGrad,depth(1:end-1),'color',clrs(iSpecies,:),'LineStyle','--')
                else if mod(iSpecies,4)==3
                        plot(relGrad,depth(1:end-1),'color',clrs(iSpecies,:),'LineStyle','-.')
                    else
                        plot(relGrad,depth(1:end-1),'color',clrs(iSpecies,:),'LineStyle',':')
                    end
                end
            end
        end
        
        set(gca,'YDir','reverse');                                                % reverse the y-axis, since depth should increase from top to bottom
        set(gcf,'Position',get(0,'ScreenSize'));                                  % maximize the figure
        set(gcf,'PaperPosition',[3,5,24,11]);                                     % set the dimensions on paper when the plots are printed
        set(gca,'FontSize',18);                                                    % set the font size
        legend(speciesName,'NumColumns',2,'Location','Best','FontSize',12)
        xlabel('\Delta C/C','FontSize',16)
        ylabel('Depth [cm]','FontSize',16)
        title(sprintf('Relative Concentration Gradient (t=%f)',time(iTimes)),'FontSize',16,'FontWeight','bold')
    end
end

end

