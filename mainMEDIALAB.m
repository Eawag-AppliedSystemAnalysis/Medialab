%**************************************************************************
%**************************************************************************
%           MEDIALAB (Modeling Early DIAgenesis using MATLAB) v.(1.1)
%                               February 2022
% MEDIALAB is an early diagenesis model, which calculates concentrations,
% fluxes of chemical species and rates of all the biogeochemical pathways 
% at each depth of aquatic sediments for a specific time period. 
% A system of partial differential equations corresponding to early
% diagenesis equations are automatically generated through MATLAB’s symbolic
% programming capabilities and solved using MATLAB’s built-in solver pdepe 
% to evaluate temporal and spatial distribution of chemical species.
% MEDIALAB is executed through MATLAB home screen through following command
%
%     [simValues, depth, time] = mainMEDIALAB
%
% Execution of MEDIALAB requires an active installation of MATLAB 
% (Version 7.6 release R2008a or later ). It is recommended to allow at 
% least 1GB of space on the hard drive for the model output and 1GB of 
% contiguous random-access memory for the initialization routine.
%**************************************************************************
%**************************************************************************


function [simValues, depth, time] = mainMEDIALAB
% mainMEDIALAB is the core script file of MEDIALAB which is executed from 
% MATLAB’s home screen. Home directory must be the same address as the 
% folder containing the contents of MEDIALAB files unless it has been saved 
% in the similar directory. It is in this script that userMEDIALAB.m,
% autoMEDIALAB and pdepe solvers are called and the final solution, sol, is
% saved in the resultMEDIALAB.mat which involves all the simulated 
% concentrations, depth and time values in a MATLAB matrix format. After 
% termination of simulation, variables of depth, time and simValues are 
% saved in the workspace. They are used for plotting and reaction rate 
% calculations through postprocessMEDIALAB.m.

% This function solves PDEs of early diageneis in lake sediments (here,
% Lake Baldegg in Switzerland)
% It uses 'pdepe' solver of MATLAB which solves 1-D elliptic PDEs
%
%
% output arguments:
%   simValues          - stores the simulated concentrations of all the 
%                        species at all depths and times. Its size is
%                        {nSpecies*1} cells where each cell stores
%                        concentrations of the species at nX depths and nT 
%                        times i.e. [nT*nX] metrices.
%   depth              - the spatial interval in [cm] with nX descritized 
%                        depths
%   time               - the temporal interval in [yr] with nT descritized
%                        times

tic;                                                                       % starts stopwatch timer                                                               % clears all the varibles in the workspace memory
close all;                                                                 % closes all the matlab figures
addpath('src','inputs')

% calling function userMEDIALAB which has been provided by the user and
% contains spatial and temporal domains and reaction network
%
% output arguments:
%   x                   - spatial interval i.e. [0,L] with nX
%                         descritized depths
%   t                   - temporal domain i.e. [0,T] with nT decritized
%                         times
%   speciesName         - cell array containing strings of names of
%                         species. Name of species end with (s) or (aq)
%                         as a flag variable to distinuish between solid or
%                         aqeouse species respectively.
%   inputFile           - string variable of name of the input file. It is
%                         a text file and has extention of .txt
%   ICFile              - string variable of name of the MATLAB file
%                         containing the initial conditions. '' if the
%                         initial conditions shouls all be set to 0
%   r                   - cell array with elements of reactions kinetic
%                         rates expresses in a string format. It has the
%                         sieze of {numRxns*1} where numRxns is number of 
%                         reactions
%   R                   - cell vector with elements of symbolic terms of
%                         reactions. Its size is {numRxns*1} 
%   isSolidReaction - vector containg 0 or 1 values when reactions have 
%                         solute or solid species units.
[x, t, speciesName, inputFile, ICFile r, R, isSolidReaction] = userMEDIALAB;  

% calling function autoMEDIALAB which automatically calculates components
% of pdepe solver i.e. vector functions of transport, reaction and boundary 
% conditions using symbolic programming capability of MATLAB
%
% output arguments:
%   isSolidSpecies      - vector of species phases i.e. 1 for solids 
%                         and 0 for solutes. It has the size of nSpecies
%   stoichiometrix      - matrix of stoichimetric coefficients. It has
%                         the size of [nSpecies*nReactions] where numRxns is
%                         the number of reactions in the reaction network
%   advection           - vector of function handle for the advection
%                         of species. It has the size of [nSpecies*1]
%   diffusion           - vector of function handle for the diffusion
%                         of species. It has the size of [nSpecies*1] 
%   reaction            - vector of function handle for reaction rates.
%                         It has the size of [1*nReactions]
%   plBC                - vector of function handle of the upper
%                         boundary condition.It has the size of 
%                         [nSpecies*1]
%   prBC                - vector of function handle of the lower
%                         boundary condition.It has the size of 
%                         [nSpecies*1]
[isSolidSpecies, stoichiometrix, advection, diffusion, reaction, plBC, prBC, parName, parValue] = autoMEDIALAB(speciesName, inputFile, r, R, isSolidReaction);


nSpecies = length(speciesName);    
nReactions = length(R);
depth = x;
time  = t;



%display starting time of solving PDE
rt_time = clock;                                    
fprintf('\nThe starting time is %.0f/%.0f/%.0f %.0f:%.0f:%.0f\n', rt_time(1),...
        rt_time(2), rt_time(3), rt_time(4), rt_time(5), rt_time(6));

    
% open logfile for storage of data at timestep level    
fid = fopen('outputs/compLevelOutput.txt', 'w'); % open file in overwrite mode
fprintf(fid,'0'); % overwrite previous data
fclose(fid);
fid = fopen('outputs/compLevelOutput.txt', 'a'); % reopen the file in append mode
disp('solving PDE ');


% calls the pdepe solver and stores the result in 'sol' which is
% [nT*nX*nSpecies] metrice. 3 functions of pdeFun, icFun and bcFun are
% called automatically by pdepe. 
m = 0;         % m value depends on the symmetry of the problem and for 1D early diagenesis problem it equals to zeo 
% options = odeset('AbsTol',1E-200,'RelTol',1E-2,'InitialStep',1E-4, 'NonNegative',ones(size(speciesName)),'BDF','on','MaxOrder',1);
% jacobian = eye(length(isSolidSpecies));
% options = odeset('AbsTol',1E-200,'RelTol',1E-6,'InitialStep',1E-4,'JPattern',jacobian);
w = str2func(strcat('@(t)',parValue{ismember(parName,'w')}));
maxStep = 1*min(diff(depth))/w(time(1))*5*100 
options = odeset('AbsTol',1E-200,'RelTol',1E-2,'InitialStep',1E-8,'MaxStep',maxStep);
sol = pdepe(m,@pdeFun,@icFun,@bcFun,x,t,options);    


%display ending time of simulation
toc;                                     
rt_time = clock;                                  
fprintf('The ending time is %.0f/%.0f/%.0f %.0f:%.0f:%.0f\n', rt_time(1),...
         rt_time(2), rt_time(3), rt_time(4), rt_time(5), rt_time(6));
     
     
% Close the logfile
fclose(fid);


% Save results and check for negative concentrations
simValues = cell(nSpecies, 1);  %creates an nSpecies*1 cell array of empty matrices for the concentration
hasNegValues = 0;
for iSpecies = 1:nSpecies
    simValues{iSpecies} = sol(:,:,iSpecies);   
    if min(min(simValues{iSpecies}))<0
        hasNegValues = 1;
        fprintf('Species %s has negative concentrations\n', speciesName{iSpecies})
    end
end      
if ~hasNegValues
    fprintf('No negative concentrations\n')
end                                           
save('outputs/resultMEDIALAB.mat','simValues','sol','time','depth');                              %saves simulation results in resultMEDIALAB matrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdeFun constructs c,f,s that are column vectors of time gradient
% coefficients, transport and reaction
%********************* TRANSPORT AND REACTION ****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,f,s] = pdeFun(x,t,u,DuDx)
    
    % always print time for a certain depth to log file
    if abs(x-depth(2)/2)<0.01
        fprintf('t=%i\n',t)
    end
    
    % coupling of the partial derivatives with respect to time is 
    % restricted to multiplication by a diagonal matrix c (x,t,u,?u/?x). 
    % In no-adsorption diagenesis problem, it equals to one
    c = ones(nSpecies,1);
    
    valAdvection = zeros(1,nSpecies);
    valDiffusion = zeros(1,nSpecies);
    %Transport vector function
    for jSpecies=1:nSpecies
        valAdvection (jSpecies) = advection{jSpecies}(u,t);
        valDiffusion (jSpecies) = diffusion{jSpecies}(DuDx,x,t);
    end
    f = valAdvection + valDiffusion;
    f=f';
    
    
    %Reaction vector function
    valReaction = zeros(1,nReactions);
    for jReaction=1:nReactions
        valReaction(jReaction) = reaction{jReaction}(u);
    end
    s = real(stoichiometrix*(valReaction'));
    
    
    % Print a line to the logfile if necessary
    [maxRatio indMax] = max(s'./u');
    if maxRatio>100000 % only store the max if it is larger than 1000
        fprintf(fid,'%e\t%f\t%i\t%e\n',t,x,indMax,maxRatio);
    else % if it is not larger, just store time for one of the depths
        if abs(x-depth(2)/2)<0.01
            fprintf(fid,'%e\t%f\t%i\t%e\n',t,x,0,0);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculates initial condition vector function for pdepe solver
%************************** INITIAL CONDITION ****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u0 = icFun(x)                            
    u0 = zeros(nSpecies,1)+1/(x+1)*1E-10; %setting all the values of the IC exactly to 0 or a fixed small number results in an error from the pdepe solver   
    if ~isempty(ICFile)
        load(ICFile, 'icValues');
        load(ICFile, 'xic');
        for jSpecies = 1:length(icValues) % Do not use nSpecies, to avoid an error in case a new species is being added
            u0(jSpecies) = interp1(xic, icValues{jSpecies}, x);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bcFun calculates boundary condition vector function at lower and upper
% boundaries
%************************* BOUNDARY CONDITION ****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pl,ql,pr,qr] = bcFun(xl,ul,xr,ur,t)   

    %Upper bounday 
    for jSpecies=1:nSpecies
        pl(jSpecies) = plBC{jSpecies}(ul,t);
    end
%     ql = isSolidSpecies;
    ql = ones(size(isSolidSpecies));
    ql=ql'; 
    pl=pl'; 
    
    %Lower bounday
    for jSpecies=1:nSpecies
        pr(jSpecies) = prBC{jSpecies}(ur,t);
    end
    pr=pr'; 
    qr = ones(nSpecies,1);
    %*************************************************************************%
end


end
