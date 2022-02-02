% This file can be used to visualize the output of mainMEDIALAB. It
% generates a number of plots of concentrations, reaction rates, mass
% balance etc. as a function of time and/or depth. Options regarding these
% figures can be set in plotOptionsMEDIALAB.
%
% DEPENDENCIES:
% - reads plotOptionsMEDIALAB to determine what to plot
% - runs userMEDIALAB and autoMEDIALAB, which are used to calculate
%   reaction rates and BC conditions (fluxes). If the userMEDIALAB file
%   does not correspond to the simValues which are being visualized, the
%   reaction rates and mass balances will be wrong.

function postprocessMEDIALAB(resultfile)
%
% input arguments:
%   simValues          - stores the simulated concentrations of all the 
%                        species at all depths and times. Its size is
%                        {nSpecies*1} cells where each cell stores
%                        concentrations of the species at nX depths and nT 
%                        times i.e. [nT*nX] metrices.
%   depth              - the spatial interval in [cm] with nX descritized 
%                        depths
%   time               - the temporal interval in [yr] with nT descritized
%                        times

close all;
addpath('src', 'inputs', 'outputs')
if nargin == 0
    load ('resultMEDIALAB')
else
    load (resultfile)
end
[pOptions] = plotOptionsMEDIALAB;
[~, ~, speciesName, inputFile, ICFile, r, R, isSolidReaction, rNames, massBalanceElements] = userMEDIALAB;
[isSolidSpecies, stoichiometrix, advection, diffusion, reaction, plBC, prBC, parName, parValue] = autoMEDIALAB(speciesName, inputFile, r, R, isSolidReaction);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************SETTING UP SOME VARIABLES****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum indices
% ------------------------------------------------------------------------
nSpecies = length(speciesName);    
nReactions = length(R);
nTimes = min(length(time), size(simValues{1},1)); % To work even if the solver aborted before finishing
nDepths = length(depth);

% Variables from inputMEDIALAB
% ------------------------------------------------------------------------
porosity = str2num(parValue{ismember(parName,'porosity')});
rho = str2num(parValue{ismember(parName,'rho')});
w = str2func(strcat('@(t)',parValue{ismember(parName,'w')}));
kEqIronSulfide = str2num(parValue{ismember(parName,'kEqIronSulfide')});
ph = str2num(parValue{ismember(parName,'ph')});
kEqViv = str2num(parValue{ismember(parName,'kEqViv')});
kEqMnCarbonate = str2num(parValue{ismember(parName,'kEqMnCarbonate')});
kEqHCO3CO3 = str2num(parValue{ismember(parName,'kEqHCO3CO3')});
seasonality = str2num(parValue{ismember(parName,'seasonality')});
convFactor = rho*(1-porosity)/porosity;

% Calculate reaction rate values
% ------------------------------------------------------------------------
if (pOptions.doPlot.reactions||pOptions.doPlot.fractionFigures) % Only do if needed (is quite time-consuming)
    ratesValues = cell(nReactions,1);
    for iReaction=1:nReactions
        ratesValues{iReaction} = zeros(nTimes,nDepths);
    end
    for iTime=1:nTimes
        for iDepth=1:nDepths
            for iSpecies=1:nSpecies
                u(iSpecies)=simValues{iSpecies}(iTime,iDepth);
            end
            for iReaction=1:nReactions
                ratesValues{iReaction}(iTime,iDepth) = reaction{iReaction}(u);   % time and space-dependent calculation of reaction rates
            end
        end
    end
    if pOptions.doPlot.fractionFigures
        for iReaction=1:nReactions
            for iTime=1:nTimes
                rateValuesTemporal(iReaction,iTime) = trapz(depth,ratesValues{iReaction}(iTime,:));        % reaction rates are integrated with respect to depth to obtain temporal rate values
            end
        end
    end
end

% Calculate fluxes
% ------------------------------------------------------------------------
if pOptions.doPlot.fractionFigures
    % Extract concentrations at upper/lower boundary
    for iSpecies=1:nSpecies
        ul(iSpecies,:) = simValues{iSpecies}(:,1);
        ur(iSpecies,:) = simValues{iSpecies}(:,nDepths);
    end
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
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************IDENTIFYING AND PLOTTING STEEPEST GRADIENTS******************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.relativeConcentrationGradients
    plotGradients(simValues, depth, time, [1, round(nTimes/2), nTimes], speciesName)
    plotGradients(simValues, depth, time, [], speciesName)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************PLOTTING EVOLUTION OF TIMESTEP*************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.timesteps
    try %might not work if simValues, depth and time are from an old simulation
        fileSize = dir('compLevelOutput.txt');
        fileSize = fileSize.bytes;
        if fileSize < 50000000 % too large files overwhelm MATLAB (memory)
            times = load('compLevelOutput.txt');
            figure;
            title('Evolution of timestep','FontSize',16,'FontWeight','bold');
            hold on
            [uTimes, nonDuplInd] = unique(times(:,1));
            plot((uTimes(1:(end-1))+uTimes(2:end))/2,log(diff(uTimes))/log(10),'b')
            plot(times(:,1),log(1./times(:,4))/log(10),'.r')
            maxStep = 1*min(diff(depth))/w(times)
            plot(times(:,1),log(maxStep)/log(10),'k')
            legend('timestep','1/max(s/u)','Max allowed timestep','Location','Best');
            xlabel('time [yr]','FontSize',16);
            ylabel('log(dt)','FontSize',16);
            
            figure;
            plot(times(:,1),times(:,3),'.k')
            ylim([0 nSpecies])
            xlabel('time [yr]','FontSize',16);
            ylabel('species number','FontSize',16);
            title('timestep-limiting species','FontSize',16,'FontWeight','bold');
        end
    catch
        % Do nothing
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%********************PLOTTING MEASURED AND MODELED PROFILES***************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.measuredProfiles
    
    % Dissolved species
    % --------------------------------------------------------------------
    figure;
    subplot(2,4,1);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'o2(aq)', {'o2(aq)'});    
    subplot(2,4,2);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'no3(aq)', {'no3(aq)'});
    subplot(2,4,3);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'so4(aq)', {'so4(aq)'});
    subplot(2,4,4);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'ch4(aq)', {'ch4(aq)'});
    subplot(2,4,5);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'nh4(aq)', {'nh4(aq)'});
    subplot(2,4,6);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'mn(aq)', {'mn(aq)'});
    subplot(2,4,7);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'fe(aq)', {'fe(aq)'});
    subplot(2,4,8);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'hpo4(aq)', {'hpo4(aq)'});
    
    
    % PLOTTING OM POOLS AND TOC PROFILES
    % --------------------------------------------------------------------
    figure;
    subplot(1,3,1);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TOC(s)', {'om1(s)'});
    subplot(1,3,2);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TOC(s)', {'om2(s)'});
    subplot(1,3,3);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TOC(s)', {'om3(s)'});
    figure;
    plotMEDIALAB(speciesName,simValues,depth,'fieldData.xlsx','TOC(s)', {'om1(s)','om2(s)','om3(s)'});
    
    
    % PLOTTING MnO2 POOLS AND TOTAL Mn(s) PROFILES
    % --------------------------------------------------------------------
    figure;
    subplot(1,4,1);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TMn(s)', {'mo_1(s)'});
    subplot(1,4,2);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TMn(s)', {'mo_2(s)'});
    subplot(1,4,3);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TMn(s)', {'mo_3(s)'});
    subplot(1,4,4);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TMn(s)', {'mo_2(s)','mo_1(s)','mo_3(s)'});
    
    
    % PLOTTING Fe(OH)3 POOLS AND TOTAL Fe(s) PROFILES
    % --------------------------------------------------------------------
    figure;
    subplot(2,3,1);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TFe(s)', {'foh_1(s)'});
    subplot(2,3,2);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TFe(s)', {'fs(s)'});
    subplot(2,3,3);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TFe(s)', {'pyrite(s)'});
    subplot(2,3,4);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TFe(s)', {'vivianite(s)'});
    subplot(2,3,5);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TFe(s)', {'foh_2(s)'});
    subplot(2,3,6);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TFe(s)', {'vivianite(s)','foh_1(s)','fs(s)','pyrite(s)','foh_2(s)'});
    
    
    % PLOTTING TOTAL P(s) PROFILE
    % --------------------------------------------------------------------
    figure;
    subplot(2,3,1);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TP(s)', {'om1(s)'});
    subplot(2,3,2);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TP(s)', {'om2(s)'});
    subplot(2,3,3);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TP(s)', {'om3(s)'});
    subplot(2,3,4);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TP(s)', {'vivianite(s)'});
    subplot(2,3,5);
    plotMEDIALAB(speciesName,simValues,depth, 'fieldData.xlsx', 'TP(s)', {'vivianite(s)','om1(s)','om2(s)','om3(s)'});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************PLOTTING TEMPORAL EVOLUTION OF REACTION PROFILES*************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.reactions
    for iReaction = 1:nReactions
        if isSolidReaction(iReaction)
            reactionName{iReaction} = sprintf('%i (s): %s',iReaction,rNames{iReaction});
        else
            reactionName{iReaction} = sprintf('%i (aq): %s',iReaction,rNames{iReaction});
        end
    end
    useLogscale = false;
    plotAgainstDepth('Reaction Rates', reactionName, ratesValues, depth, {'[\mumol/g/yr]', '[\mumol/cm^3/yr]'}, useLogscale)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**********PLOTTING TEMPORAL EVOLUTION OF CONCENTRATION PROFILES**********%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.concentrations
    useLogscale = true;
    plotAgainstDepth('Concentrations', speciesName, simValues, depth, {'\mumol/g', '\mumol/cm^3'}, useLogscale)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*******PLOTTING TEMPORAL EVOLUTION OF SATURATION INDICES PROFILES********%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.saturationIndices
    useLogscale = true;
    speciesNameSat = {'Iron Sulfide','Vivianite','Mn Carbonate'};
    simValuesSat{1} = simValues{8}.*simValues{9}/(kEqIronSulfide*ph^1);
    simValuesSat{2} = simValues{8}.^3.*simValues{14}.^2/kEqViv;
    simValuesSat{3} = simValues{7}.*simValues{19}*kEqHCO3CO3/(ph*kEqMnCarbonate);
    plotAgainstDepth('Saturation Indices', speciesNameSat, simValuesSat, depth, {'-', '-'}, useLogscale)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************* GENERATING ALL FRACTION PLOTS ***********************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.fractionFigures
    
    % Initializationn Part 1
    % --------------------------------------------------------------------
    X=time(1:nTimes)';
    Z = 0;
    
    % Plotting origin fractions over time (depth-integrated)
    % --------------------------------------------------------------------
    if strcmp(pOptions.fractionProductionSelection,'all')
        pOptions.fractionProductionSelection = speciesName;
    end
    for index=1:numel(pOptions.fractionProductionSelection)
        iSpecies = find(strcmp(speciesName,pOptions.fractionProductionSelection{index}) ~= 0);
        Y=[];
        legends = {};
        for iRxn = 1:nReactions
            rTemp = max(0,stoichiometrix(iSpecies,iRxn)*rateValuesTemporal(iRxn,:));
            if isSolidReaction(iRxn)
                rTemp = rTemp*rho*(1-porosity);
            else
                rTemp = rTemp*porosity;
            end
            if min(rTemp)>0
                Y = [Y; rTemp];
                legends = {legends{:}, rNames{iRxn}};
            end
        end
        fluxInTop = max(0,fluxl{iSpecies});
        Y = [Y; fluxInTop];
        legends = {legends{:}, 'fTop'};
        
        plotFractions(sprintf('%s sources',speciesName{iSpecies}), X, Y, Z, 'time', legends, pOptions.fractionRelative)
    end


    % Plotting consumption/degradation fractions over time (depth-integrated, including OM degradation)
    % --------------------------------------------------------------------
    if strcmp(pOptions.fractionConsumptionSelection,'all')
        pOptions.fractionConsumptionSelection = speciesName;
    end
    for index=1:numel(pOptions.fractionConsumptionSelection)
        iSpecies = find(strcmp(speciesName,pOptions.fractionConsumptionSelection{index}) ~= 0);
        Y=[];
        legends = {};
        for iRxn = 1:nReactions
            rTemp = max(0,-stoichiometrix(iSpecies,iRxn)*rateValuesTemporal(iRxn,:));
            if isSolidReaction(iRxn)
                rTemp = rTemp*rho*(1-porosity);
            else
                rTemp = rTemp*porosity;
            end
            if min(rTemp)>0
                Y = [Y; rTemp];
                legends = {legends{:}, rNames{iRxn}};
            end
        end
        fluxOutTop = max(0,-fluxl{iSpecies});
        fluxOutBottom = max(0,fluxr{iSpecies});
        Y = [Y; fluxOutTop; fluxOutBottom];
        legends = {legends{:}, 'fTop', 'fBottom'};
        
        plotFractions(sprintf('%s sinks',speciesName{iSpecies}), X, Y, Z, 'time', legends, pOptions.fractionRelative)
    end

    
    % Initializationn Part 2
    % --------------------------------------------------------------------
    X=depth';
    ind = [1, round(nTimes/3), round(nTimes/3*2), nTimes]; % Determine which time steps to plot
    Z = time(ind);
    
    
    % Plotting production fractions as a function of depth at
    % four different times (including OM degradation)
    % --------------------------------------------------------------------
    if strcmp(pOptions.fractionProductionDepthSelection,'all')
        pOptions.fractionProductionDepthSelection = speciesName;
    end
    for index=1:numel(pOptions.fractionProductionDepthSelection)
        iSpecies = find(strcmp(speciesName,pOptions.fractionProductionDepthSelection{index}) ~= 0);
        Y=[];
        legendSelection = zeros(nReactions,1);
        legends = {};
        for kTime = 1:numel(ind)
            jCount = 0;
            for jRxn = 1:nReactions
                rTemp = max(0,stoichiometrix(iSpecies,jRxn)*ratesValues{jRxn}(ind(kTime),:));
                if isSolidReaction(jRxn)
                    rTemp = rTemp*rho*(1-porosity);
                else
                    rTemp = rTemp*porosity;
                end
                if max(rTemp)>0
                    jCount = jCount+1;
                    Y(jCount,:,kTime) = rTemp;
                    legendSelection(jRxn)=1;
                end
            end
        end
        for jRxn=1:nReactions
            if legendSelection(jRxn)
                legends = {legends{:}, rNames{jRxn}};
            end
        end
        plotFractions(sprintf('%s sources',speciesName{iSpecies}), X, Y, Z, 'depth', legends, pOptions.fractionRelative)
    end
    
    
    % Plotting consumption/degradation fractions as a function of depth at
    % four different times (including OM degradation)
    % --------------------------------------------------------------------
    if strcmp(pOptions.fractionConsumptionDepthSelection,'all')
        pOptions.fractionConsumptionDepthSelection = speciesName;
    end
    for index=1:numel(pOptions.fractionConsumptionDepthSelection)
        iSpecies = find(strcmp(speciesName,pOptions.fractionConsumptionDepthSelection{index}) ~= 0);
        Y=[];
        legendSelection = zeros(nReactions,1);
        legends = {};
        for kTime = 1:numel(ind)
            jCount = 0;
            for jRxn = 1:nReactions
                rTemp = max(0,-stoichiometrix(iSpecies,jRxn)*ratesValues{jRxn}(ind(kTime),:));
                if isSolidReaction(jRxn)
                    rTemp = rTemp*rho*(1-porosity);
                else
                    rTemp = rTemp*porosity;
                end
                if max(rTemp)>0
                    jCount = jCount+1;
                    Y(jCount,:,kTime) = rTemp;
                    legendSelection(jRxn)=1;
                end
            end
        end
        for jRxn=1:nReactions
            if legendSelection(jRxn)
                legends = {legends{:}, rNames{jRxn}};
            end
        end
        plotFractions(sprintf('%s sinks',speciesName{iSpecies}), X, Y, Z, 'depth', legends, pOptions.fractionRelative)
    end
    
    
    % Plotting fraction different species are contributing to elemental
    % content as a function of depth at four different times
    % --------------------------------------------------------------------
    if strcmp(pOptions.fractionElementalSelection,'all')
        fTemp = fields(massBalanceElements);
        pOptions.fractionElementalSelection = fTemp(1:(end-1)); % do not include Fred
    end
    for iElement=1:numel(pOptions.fractionElementalSelection)
        Y = [];
        legendEntries = {};
        jCount = 0;
        elementContent = massBalanceElements.(pOptions.fractionElementalSelection{iElement});
        for jSpecies = 1:nSpecies
            if (elementContent(jSpecies)>0)&&(isSolidSpecies(jSpecies)||pOptions.fractionElementalIncludeSolutes)
                jCount = jCount+1;
                for kTime = 1:numel(ind)
                    fTemp = max(0,elementContent(jSpecies)*simValues{jSpecies}(ind(kTime),:));
                    if isSolidSpecies(jSpecies)
                        fTemp = fTemp*rho*(1-porosity);
                    else
                        fTemp = fTemp*porosity;
                    end
                    Y(jCount,:,kTime) = fTemp;
                end
                legendEntries = {legendEntries{:}, speciesName{jSpecies}};
            end
        end
        plotFractions(pOptions.fractionElementalSelection{iElement}, X, Y, Z, 'depth', legendEntries, pOptions.fractionRelative)
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%******PLOTTING TEMPORAL CONCENTRATIONS NEAR SEDIMENT-WATER INTERFACE*****%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.concentrationsTemporal
    iDepthToPlot = [1,11,51];
    useLogscale = false;
    plotAgainstTime('Concentrations', speciesName, simValues, time(1:nTimes), depth, iDepthToPlot, {'\mumol/g', '\mumol/cm^3'}, useLogscale)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%********PLOTTING TEMPORAL FLUXES ACROSS SEDIMENT-WATER INTERFACE*********%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.fluxesTemporal
    for iSpecies=1:nSpecies
        DuDx(iSpecies,:) = (simValues{iSpecies}(:,1)-simValues{iSpecies}(:,2))/(depth(1)-depth(2));
    end
    for iSpecies=1:nSpecies
        for iTime=1:nTimes
            flux{iSpecies}(iTime,1) = diffusion{iSpecies}(DuDx(:,iTime),(depth(1)+depth(2))/2,time(iTime))*(~isSolidSpecies(iSpecies));
        end
    end
    
    useLogscale = false;
    plotAgainstTime('Solute fluxes at upper BC',speciesName, flux, time(1:nTimes), depth(2)/2, 1,{'\mumol/cm^2/yr'}, useLogscale)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%**********PLOTTING MASS BALANCE FOR MULTIPLE SPECIES*********************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.massBalances
    massBalElNames = fields(massBalanceElements);
    for n = 1:numel(massBalElNames)-1; %Do not include Fred
        plotElementalFluxes(massBalElNames{n}, massBalanceElements.(massBalElNames{n}), speciesName, plBC, prBC, simValues, time(1:nTimes), depth, isSolidSpecies, '\mumol/cm^2/yr', porosity, rho)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%***************PLOTTING TEMPORAL CHANGE OF FRED**************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if pOptions.doPlot.fRed
    o2InputStr = strcat('@(t)',parValue{ismember(parName,'BC_o2')});
    o2InputStr = strrep(o2InputStr, 'seasonality', num2str(seasonality));
    o2Input = str2func(o2InputStr);
    plotFred(massBalanceElements.Fred, speciesName, plBC,simValues, time(1:nTimes), isSolidSpecies, 'g O_2/m^2 d', porosity, rho, o2Input)
end