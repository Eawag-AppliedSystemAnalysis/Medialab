function plotFractions(titleText,X,Y,Z,xVar,legendEntries,relativeFraction)

% plotFractions plots the relative fractions of the  Y values for all
% instances listed in Z. Values in Y are automatically scaled to give a
% total of 1. Fractions can e.g. be plotted against depth or time. If Z 
% contains more than one instance, subplots are automatically generated 
% for each instance in Z. Values in Y 
%
%
% input arguments are:
%   titleText   	- title of the figure
%   X               - Vector of values for the x-axis (e.g. depth or time)
%   Y               - 2D or 3D matrix of the components that should be
%                     plotted. Dimensions [length(X)*nComponents*nPlots].
%   Z               - Vector of z-values of each instance (e.g. time or 
%                     depth) of
%                     length nPlots
%   xVar            - string indicating what the x variable is. This 
%                     affects the text of labels and titles. Currently
%                     only 'depth' and 'time' are implemented.
%   legendEntries   - cell array of strings used for the legend. Should
%                     containt nComponents entries.
%  relativeFraction - if true, concentrations/reaction rates are scaled to 
%                     give a total of 100% at each depth/time, if false,
%                     they are not scaled (i.e. the absolute
%                     concentrations/reaction rates are plotted (stacked))
% ------------------------------------------------------------------------%


% Initialization
% ------------------------------------------------------------------------
nComponents = size(Y,1);
if ndims(Y)==3
    nPlots = size(Y,3);
else
    nPlots = 1;
end
if numel(Y)==0
    disp(sprintf('No entries for Figure: %s\nFigure not plotted',titleText))
else
    figure;
    if nPlots>1
        set(gcf,'Position',get(0,'ScreenSize'));                          %maximize the figure
        set(gcf,'PaperPosition',[3,5,24,11]);                             %set the dimensions on paper when the plots are printed
    end
    for iComp = 1:nComponents
        assert(min(min(Y(iComp,:,:)))>=0,sprintf('Plotting negative values not allowed: Component %s has negative values',legendEntries{iComp}))
    end


    
    % Scaling of Y to give a total of 1 (=100%)
    % ------------------------------------------------------------------------
    if relativeFraction
        Ytot = zeros(1,size(Y,2),nPlots);
        for iComp=1:nComponents
            Ytot = Ytot+Y(iComp,:,:);
        end
        for iComp=1:nComponents
            Y(iComp,:,:) = Y(iComp,:,:)./Ytot;
        end
    else
        yLimMax = max(max(sum(Y,1)));
    end
    
    
    % Generate the Plot
    % ------------------------------------------------------------------------
    for iPlot = 1:nPlots
        subplot(ceil(sqrt(nPlots)),ceil(sqrt(nPlots)),iPlot);
        area(X,Y(:,:,iPlot)');
        if relativeFraction
            ylim([0 1]);
        else
            ylim([0 yLimMax])
        end
        
        
        % label everything
        % --------------------------------------------------------------------
        
        % legend
        if iPlot==1
            legend(legendEntries,'Location','Best');
        end
        
        % title and xlabel + ylabel
        if nPlots==1
            if relativeFraction
                completeTitle = sprintf('relative fraction %s',titleText);
            else
                completeTitle = sprintf('absolute fraction %s',titleText);
            end
        end
        switch xVar
            case 'depth'
                xlabel('depth (cm)','FontSize',16);
                if nPlots>1
                    if relativeFraction
                        completeTitle = sprintf('relative fraction %s (%0.2f years)',titleText,Z(iPlot));
                    else
                        completeTitle = sprintf('absolute fraction %s (%0.2f years)',titleText,Z(iPlot));
                    end
                end
                if relativeFraction
                    ylabel('Fraction');
                else
                    ylabel('\mumol/cm^3 (bulk)/yr');
                end
                    
            case 'time'
                xlabel('time (year)','FontSize',16);
                if nPlots>1
                    if relativeFraction
                        completeTitle = sprintf('relative fraction %s (%0.1f cm)',titleText,Z(iPlot));
                    else
                        completeTitle = sprintf('absolute fraction %s (%0.1f cm)',titleText,Z(iPlot));
                    end
                end
                if relativeFraction
                    ylabel('Fraction');
                else
                    ylabel('\mumol/cm^2/yr');
                end
                
            otherwise
                % no xlabel
                if nPlots>1
                    if relativeFraction
                        completeTitle = sprintf('relative fraction %s (%0.1f)',titleText,Z(iPlot));
                    else
                        completeTitle = sprintf('absolute fraction %s (%0.1f)',titleText,Z(iPlot));
                    end
                end
        end
        title(completeTitle,'FontSize',16,'FontWeight','bold');
        
        
    end
    
end

