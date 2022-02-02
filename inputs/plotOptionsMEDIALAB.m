% plotOptionsMEDIALAB allows the user to select which plots should be
% generated and to set some options concerning the plots


function  [pOptions] = plotOptionsMEDIALAB  
% 
% Output arguments:
%                   -pOptions is a struct containing various settings used
%                   for plotting


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************** SELECTION OF PLOT TYPES ****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots relative concentration gradients (dC/dz)/C(z) as a function of
% depth (at initial time, half the simulation time, and at the end). These
% can e.g. be used to diagnose whether the difference in resolution between
% top and bottom of the modeled sediment slab is appropriate. Since the
% gradients are calculated per cell, rather than per cm, the highest
% gradients observed in the upper part of the model should be similar to
% those in the lower part.
pOptions.doPlot.relativeConcentrationGradients = false;

% Plots how the size of the timestep changes over time
pOptions.doPlot.timesteps = false;

% Plots measured profiles of solute and solid species from fieldData.xlsx
pOptions.doPlot.measuredProfiles = false;

% Plots modeled reaction rates for several times as a function of depth
pOptions.doPlot.reactions = true;

% Plots modeled concentrations for several times as a function of depth
pOptions.doPlot.concentrations = true;

% Plots Saturation Indices for several times as a function of depth
pOptions.doPlot.saturationIndices = false;

% Plots figures showing which fractions of certain elements or species are 
% contained in which species or consumed/produced by which reaction
pOptions.doPlot.fractionFigures = false;

% Plots concentrations as a function of time for three (shallow) depths
pOptions.doPlot.concentrationsTemporal = false;

% Plots fluxes for each species as a function of time
pOptions.doPlot.fluxesTemporal = false;

% Plots mass balance for several elements over time. This can also be used
% to check whether changes made to the stoichiometry of some reactions are
% correct. If the mass balances do not hold, probably the stoichiometry of
% a reaction is wrong.
pOptions.doPlot.massBalances = false;

% Plots F_red over time
pOptions.doPlot.fRed = false;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************** SETTING FURTHER OPTIONS ****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determines whether fractions are scaled to always give 100% or not
pOptions.fractionRelative = true;

pOptions.fractionProductionSelection = {'no3(aq)','fe(aq)','s2(aq)','fs(s)','vivianite(s)'};
% pOptions.fractionProductionSelection = {'all'};

pOptions.fractionProductionDepthSelection = {'nh4(aq)','fe(aq)','s2(aq)','fs(s)','vivianite(s)'};
% pOptions.fractionConsumptionDepthSelection = {'all'};


pOptions.fractionConsumptionSelection = {'om1(s)','o2(aq)','no3(aq)','fe(aq)','s2(aq)','fs(s)','vivianite(s)'};
% pOptions.fractionConsumptionSelection = {'all'};


pOptions.fractionConsumptionDepthSelection = {'om1(s)','fe(aq)','nh4(aq)','s2(aq)','fs(s)','vivianite(s)'};
% pOptions.fractionConsumptionDepthSelection = {'all'};


pOptions.fractionElementalSelection = {'Fe','S'};
% pOptions.fractionElementalSelection = {'all'};
pOptions.fractionElementalIncludeSolutes = false;
