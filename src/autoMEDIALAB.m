% Transport, reaction, and boundary condition vectors 
% are automatically generated via autoMEDIALAB.m in form of MATLAB 
% function handles and are passed to mainMEDIALAB.m. 

function  [isSolidSpecies, stoichiometrix, advection, diffusion, reaction, plBC, prBC, parName, parValue] = autoMEDIALAB(speciesName,inputFile, r, R, isSolidReaction)

% This function :
%   1) Reads the model parameters from inputMEDIALAB.txt
%   2) constructs reaction network stoichiometrics matrix as well as 
%      reaction transport and upper boundary condition components vector 
%      functions in the form of MATLAB function handles.
%
% input arguments:
%   speciesName         - cell array containing strings of names of
%                         species. Name of species end with (s) or (aq)
%                         as a flag variable to distinuish between solid or
%                         aqeouse species respectively.
%   inputFile           - string variable of name of the input file. It is
%                         a text file and has extention of .txt
%   r                   - cell array with elements of reactions kinetic
%                         rates expresses in a string format. It has the
%                         sieze of {numRxns*1} where numRxns is number of 
%                         reactions
%   R                   - cell vector with elements of symbolic terms of
%                         reactions. Its size is {numRxns*1} 
%   isSolidReaction     - vector containg 0 or 1 values when reactions have 
%                         solute or solid species units.
%
%
% output arguments:
%   isSolidSpecies       - vector of species phases i.e. 1 for solids 
%                         and 0 for solutes. It has the size of nSpecies
%   stoichiometrix       - matrix of stoichimetric coefficients. It has
%                         the size of [nSpecies*nReactions] where numRxns is
%                         the number of reactions in the reaction network
%   advection           - vector of function handle for the advection term
%                         of species. It has the size of [nSpecies*1]
%   diffusion           - vector of function handle for the diffusion
%                         of species. It has the size of [nSpecies*1] 
%   reaction            - vector of function handle for reaction rates.
%                         It has the size of [nReactions*1]
%   plBC                - vector of function handle of the upper
%                         boundary condition.It has the size of 
%                         [nSpecies*1]
%   prBC                - vector of function handle of the lower
%                         boundary condition.It has the size of 
%                         [nSpecies*1]
%   parName             - vector of names of parameters specified in inputMEDIALAB
%   parValue            - vector of values of parameters specified in inputMEDIALAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%******************************1st BLOCK**********************************%
%************************** READING INPUT DATA ***************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('**********************');
disp('READING THE INPUT FILE');
fid            = fopen(inputFile,'r');                                    % reading the parameters from input file
C              = textscan(fid,'%s%s%s%s');                                    % %s ist ein Formatierungs Befehl
rateConsNames  = C{1};                                                    % extracting the names of parameters
rateConsVals   = C{2};                                                    % extracting the values of parameters
transRateConsVals = C{3};
seasonality    = C{4};
fclose(fid);                                                              % close file


% Construct the parameter strings for those parameters which change over
% time. This is indicated in the inputMEDIALAB: if after the name of any
% parameter, only a single number is given, it is assumed to be constant
% during all of the simulation time. If the number is followed by an array
% of years and values, the parameter changes over time. Between the given
% time/value pairs, linear interpolation is used. Before the first and
% after the last time/value pair, the parameter is assumed to be constant.
for i=1:length(rateConsVals)
    if ~isempty(transRateConsVals{i})
        temp = str2num(transRateConsVals{i});
        years = temp(:,1);
        vals = temp(:,2);
        rateConsVals{i} = strcat(rateConsVals{i},'*(');
        if vals(1)~=0 %in case this is zero, the term will be 0 anyway
            rateConsVals{i} = strcat(rateConsVals{i},'(t<',num2str(years(1)),')*',num2str(vals(1))); % before 1st date, values are constant
        end
        for j = 1:(length(years)-1)
            if (vals(j+1)-vals(j))==0 %in case this is zero, the term can be shorter
                rateConsVals{i} = strcat(rateConsVals{i},'+((t>=',num2str(years(j)),')&&(t<',num2str(years(j+1)),'))*',num2str(vals(j)));
            else
                rateConsVals{i} = strcat(rateConsVals{i},'+((t>=',num2str(years(j)),')&&(t<',num2str(years(j+1)),'))*(',num2str(vals(j)),'+',num2str((vals(j+1)-vals(j))/(years(j+1)-years(j))),'*(t-',num2str(years(j)),'))');
            end
        end
        if vals(end)~=0 %in case this is zero, the term will be 0 anyway
            rateConsVals{i} = strcat(rateConsVals{i},'+(t>=',num2str(years(end)),')*',num2str(vals(end))); % after last date, values are constant
        end
        rateConsVals{i} = strcat(rateConsVals{i},')');
    end
    rateConsVals{i}
end


for i=1:length(rateConsVals)
    disp(strcat(rateConsNames(i), ' = ',rateConsVals(i)));                % printing on the screen names and values of parameters
end

% Define some variables that are needed later on
convFactor = str2num(rateConsVals{ismember(rateConsNames,'rho')})...
    *(1-str2num(rateConsVals{ismember(rateConsNames,'porosity')}))...
    /str2num(rateConsVals{ismember(rateConsNames,'porosity')})
porosity = str2num(rateConsVals{ismember(rateConsNames,'porosity')});
dblThickness = str2num(rateConsVals{ismember(rateConsNames,'dblThickness')});
nSpecies = int16(length(speciesName));
nReactions = int16(length(R));

disp('**********************');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************************2nd BLOCK**********************************%
%********************** EXTRACTING PHASES OF SPECIES *********************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ruft species userMediaLab auf zählt sie und checkt ob species aq oder s
% und wandelt dies in 0 bzw 1 um
for iSpecies=1:length(speciesName)
    if (isempty(strfind(speciesName{iSpecies},'(aq)')));
        isSolidSpecies(iSpecies) = 1;                                               % isSolidSpecies equals to 1 if species is a solid
        speciesName{iSpecies} = strrep(speciesName{iSpecies},'(s)','');
    else
        isSolidSpecies(iSpecies)=0;                                                 % isSolidSpecies equals to 0 if species is a solute
        speciesName{iSpecies} = strrep(speciesName{iSpecies},'(aq)','');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************************3rd BLOCK**********************************%
%**************** CALCULATING REACTION NETWORK STOICHIMETRICS ************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CALCULATING REACTION NETWORK STOICHIOMETRICS');
stoichiometrix = zeros(nSpecies,length(R));
for iSpecies=1:nSpecies
    for jReaction=1:nReactions
        c=coeffs(R{jReaction},speciesName{iSpecies});                      % extracting coeffiecints of each reaction
        if (length(c)>1) % i.e. the species i is part of reaction j
            if ((isSolidSpecies(iSpecies)&&isSolidReaction(jReaction))||((~isSolidSpecies(iSpecies))&&(~isSolidReaction(jReaction))))
                stoichiometrix(iSpecies,jReaction) = -c(2);                % if both species and reaction have same units, the coefficient wouldn't be corrected by conversion factor otherwise conFactor will be imposed as follows
            elseif (isSolidSpecies(iSpecies)&&(~isSolidReaction(jReaction)))
                stoichiometrix(iSpecies,jReaction) = -c(2)/convFactor;
            else
                stoichiometrix(iSpecies,jReaction) = -c(2)*convFactor;
            end;
        else
            stoichiometrix(iSpecies,jReaction) = 0;                        % if species is not incorporated in the reaction, then the coefficient equals to zero
        end
    end
end
disp(num2str(double(stoichiometrix)));
disp('**********************');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************************4th BLOCK**********************************%
%************** CONSTRUCTING FUNCTION HANDLE OF REACTION RATES ***********%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CONSTRUCTING FUNCTION HANDLE OF REACTION RATES');
for iReaction=1:nReactions
    for jSpecies=1:nSpecies
        uIndex     = strcat('u(',strcat(num2str(jSpecies),')'));
        rnew{iReaction} = strrep(r{iReaction}, speciesName{jSpecies}, uIndex); % reaction rate terms 'r' are replaced by their numeric values u(i) where i is the index of the species in the sequential vector of diagenesis equations
        r{iReaction}    = rnew{iReaction};
    end
    for jRateCons=1:length(rateConsVals)
        rnew{iReaction} = strrep(r{iReaction}, rateConsNames{jRateCons}, char(rateConsVals(jRateCons)));
        r{iReaction}    = rnew{iReaction};
    end
    r{iReaction} = strcat('@(u)',r{iReaction});
%     r{iReaction}
    r{iReaction} = str2func(r{iReaction});
end
disp(' FUNCTION HANDLE OF REACTION RATES = '); 
r
disp('**********************');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******************************5th BLOCK**********************************%
% CONSTRUCTING FUNCTION HANDLE OF ADVECTION, DIFFUSION AND BOUNDARY CONDITIONS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('CONSTRUCTING FUNCTION HANDLE OF ADVECTION, DIFFUSION AND BOUNDARY CONDITION');

D_bio = strcat('D_b0*(1-tanh((x-H_bio)/tau_bio))*',num2str(1/(1-tanh(-str2num(rateConsVals{ismember(rateConsNames,'H_bio')})/str2num(rateConsVals{ismember(rateConsNames,'tau_bio')})))));

for iSpecies=1:nSpecies
    
    % Construct Strings which will later be converted to functions
    advString {iSpecies} = strcat('@(u,t)-w*u(',num2str(iSpecies),')');
    prString{iSpecies}   = strcat('@(ur,t)w*ur(',num2str(iSpecies),')'); %was previously multiplied additionally with the porosity
    if isSolidSpecies(iSpecies) % for solids
        plString{iSpecies}   = strcat('@(ul,t)(BC_',speciesName{iSpecies},')/',num2str(convFactor*porosity));
        difString {iSpecies} = strcat('@(DuDx,x,t)(',D_bio,'+D_bmin)*DuDx(',num2str(iSpecies),')');
    else % for solutes
        Dmol = str2num(rateConsVals{ismember(rateConsNames,strcat('Dmol_',speciesName{iSpecies}))});
        D             = num2str(Dmol/(1-log(porosity^2)));
        plString{iSpecies}   = strcat('@(ul,t)Dmol_',speciesName{iSpecies},'*(BC_',speciesName{iSpecies},'-ul(',num2str(iSpecies),'))/',num2str(dblThickness*porosity));
        difString{iSpecies} = strcat('@(DuDx,x,t)(',D_bio,'+',D,')*DuDx(',num2str(iSpecies),')');
    end
    
    % replace the rate constant names with the actual values
    for jRateCons=1:length(rateConsVals)
        advStringNew{iSpecies} = strrep(advString{iSpecies}, rateConsNames{jRateCons}, char(rateConsVals(jRateCons)));
        advString{iSpecies}    = advStringNew{iSpecies};
        difStringNew{iSpecies} = strrep(difString{iSpecies}, rateConsNames{jRateCons}, char(rateConsVals(jRateCons)));
        difString{iSpecies}    = difStringNew{iSpecies};
        plStringNew{iSpecies}  = strrep(plString{iSpecies}, rateConsNames{jRateCons}, char(rateConsVals(jRateCons)));
        plString{iSpecies}     = plStringNew{iSpecies};
        prStringNew{iSpecies}  = strrep(prString{iSpecies}, rateConsNames{jRateCons}, char(rateConsVals(jRateCons)));
        prString{iSpecies}     = prStringNew{iSpecies};
    end
    
    % Convert the strings to function handles
    advHandle{iSpecies}    = str2func(advString{iSpecies});
    advString{iSpecies}
    difString{iSpecies}
    difHandle{iSpecies}    = str2func(difString{iSpecies});
    plHandle{iSpecies}     = str2func(plString{iSpecies});
    prHandle{iSpecies}     = str2func(prString{iSpecies});
end

% Print out the function handles
disp(' FUNCTION HANDLE OF ADVECTION = ');
advHandle
disp('**********************');

disp(' FUNCTION HANDLE OF DIFFUSION = ');
difHandle
disp('**********************');

disp(' FUNCTION HANDLE OF UPPER BOUNDARY CONDITION = '); 
plHandle
disp('**********************');

disp(' FUNCTION HANDLE OF LOWER BOUNDARY CONDITION = '); 
prHandle
disp('**********************');


%*************************************************************************%
stoichiometrix = double (stoichiometrix);
reaction  = r;
advection = advHandle;
diffusion = difHandle;
plBC      = plHandle;
prBC      = prHandle;
parName  = rateConsNames;
parValue = rateConsVals;
  
 
             
             
        