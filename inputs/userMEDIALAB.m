% userMEDIALAB is the first script that has to be modified by the user to
% provide all the necessary information of spatial and temporal domains as 
% well as biogeochemical reaction network. There are some other parameters
% such as list of species names, organic matter composition coefficients 
% and input file name that have to be included in the script. Species are 
% characterized by string variables and their assigned names contain either
% �(s)� or �(aq)� representing solid or aqueous species, respectively. This
% is an approach to characterize cloase of the species in the script when 
% building the transport matrix (as molecular diffusion is excluded in case
% of solid species) and reaction matrix ( to ensure consistency between 
% reactions and species units). As described above, variable names used for 
% the chemical species are defined as symbols to set up the reaction 
% network. 

function  [x, t, speciesName, inputFile, ICFile, r, R, isSolidReaction, rName, massBalanceElements] = userMEDIALAB  

% In this function which must be edited by the user:
%  1) spatial and temporal domains of x and t are defined
%  2) chemical composition of OM i.e carbon (Cx), nitrogen(Cy) and 
%     phosphorous (Cz) content of OM is defined. In the current model there
%     are 3 pools of OM which represent labile, less labile and refractory 
%     fractions
%  3) text file, inputFile.txt, containing the input parameters is provided
%  4) MATLAB file containing the initial conditions
%  5) name of the species that diagensis equation is solved for
%  6) user defines the reactions, R, in algebraic form applying 
%     MATLAB's symbolic programming capabilitieseach and their 
%     corresponding rates, r, in a string format . For reaction defined  
%     by R{i} in this script there are three attributes: reaction rates 
%     defined by r{i} isSolidReaction(i) and rName{i}.
%
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
%                         initial conditions should all be set to 0
%   r                   - cell array with elements of reactions kinetic
%                         rates expresses in a string format. It has the
%                         sieze of {numRxns*1} where numRxns is number of 
%                         reactions
%   R                   - cell vector with elements of symbolic terms of
%                         reactions. Its size is {numRxns*1} 
%   isSolidReaction     - vector containg 0 or 1 values when reactions have 
%                         solute or solid species units.
%   rName               - cell array containing strings of names of
%                         reactions, which are used for plotting purposes
%                         only
%   massBalanceElements - Struct defining for selected elements (e.g. C, N,
%                         Fe) in which amounts they are present in each 
%                         species. This is required for the mass balance
%                         plots. Additionally, it also contains a field
%                         Fred where it is defined how each species
%                         contributes to Fred.


% %************************** USER DEFINED *********************************
% definining the spatial resolution and total sediment depth to be modeled.
% x = linspace(0, L, L/0.1) creats a domain of L/0.1 depths in the spatial 
% interval of [0,L]. The depth resolution of 0.1cm is an adequate precision 
% for most of the diagenesis problem. Units are in [cm]. Here 20 cm of
% sediment is simulated
nGridPoints = 500;
maxDepth = 50;
curvature = maxDepth/10*5*5*5; % The larger the curvature, the larger the difference in resolution between top and bottom
x = maxDepth*(curvature.^(0:1/(nGridPoints-1):1)-1)/(curvature-1); %gradually decreasing resolution
% x = linspace(0,50,1000);

% defining time span of the simulation and the temporal resolution.
% t = linspace(0, T, nT) cteats a simuation period of T years with nT
% descritized points. Units are in [yr].
t = linspace(0,20,160)+1850; %runs for 20 years with 8 outputs per year, starting in the year 1850;
% add 0.25 to start when periodic BC are at mean



% input file for the model's parameters
inputFile = 'InputMEDIALAB.txt';
% File containing the initial conditions. The file has to contain a 
% vector xic providing the spatial grid of the initial conditions and a
% cell array icValues with the concentrations at the depths specified in
% xic. If the initial conditions should all be set to 0, ICFile can be set
% to '' (an empty string).
ICFile = 'icValues.mat';
% ICFile = '';


% Open the input file to read in the cx1, cz1, etc.
 fid            = fopen(inputFile,'r');                                    % reading the parameters from input file
 C              = textscan(fid,'%s%s%s%s');                                    % %s ist ein Formatierungs Befehl
 rateConsNames  = C{1};                                                    % extracting the names of parameters
 rateConsVals   = C{2};                                                    % extracting the values of parameters
 fclose(fid);    

% carbon content of 1st OM pool 
cx1 = str2num(rateConsVals{ismember(rateConsNames,'cx1')});
% nitrogen content of 1st OM pool
cy1 = str2num(rateConsVals{ismember(rateConsNames,'cy1')});
% phosphorous content of 1st OM pool
cz1 = str2num(rateConsVals{ismember(rateConsNames,'cz1')});
% sulfur content of 1st OM pool
cs1 = str2num(rateConsVals{ismember(rateConsNames,'cs1')});
% carbon content of 2nd OM pool
cx2 = str2num(rateConsVals{ismember(rateConsNames,'cx2')});
% nitrogen content of 2nd OM pool
cy2 = str2num(rateConsVals{ismember(rateConsNames,'cy2')});
% phosphorous content of 2nd OM pool
cz2 = str2num(rateConsVals{ismember(rateConsNames,'cz2')});
% sulfur content of 2nd OM pool
cs2 = str2num(rateConsVals{ismember(rateConsNames,'cs2')});


% species names including their phase i.e. solid or solute which is denoted
% by (s) or (aq)
speciesName = {'o2(aq)',...      %1 oxygen  
               'no3(aq)',...     %2 nitrate
               'mo_1(s)',...     %3 1st manganese oxide pool (fast)
               'foh_1(s)',...    %4 iron hydroxide pool 1
               'so4(aq)',...     %5 sulfate
               'nh4(aq)',...     %6 ammonium
               'mn(aq)',...      %7 manganese(II) 
               'fe(aq)',...      %8 iron(II)
               's2(aq)',...      % sulfide (HS)
               'fs(s)',...       %10 iron sulfide (FeS)
               'pyrite(s)',...   %11 pyrite
               'om1(s)',...      %12 1st organic matter pool (fast)
               'om2(s)', ...     %13 2nd organic matter pool (slow)
               'hpo4(aq)',...    %14 phosphate
               'ch4(aq)',...     %15 methane
               'om3(s)',...      %16 3rd organic matter pool (not reactive)
               'mo_2(s)',...     %17 2nd manganese oxide (slow)
               'vivianite(s)',...%18 vivianite Fe3(PO4)2
               'hco_3(aq)',...   %19
               'co_2(aq)',...    %20
               'n2(aq)',...      %21
               'foh_2(s)',...    %22
               'mo_3(s)'};       %23
       
% all of the species included in the reactions must be included in this
% list
syms om1 om2 om3 mo_1 mo_2 foh_1 foh_2 fs pyrite vivianite...
     o2 no3 so4 nh4 mn fe s2 ch4 hpo4 hco_3 co_2 h2o n2 h h2 mo_3
        
 
% Struct defining for selected elements (e.g. C, N, Fe) in which amounts
% they are present in each species. This is required for the mass balance
% plots. Additionally, it also contains a field Fred where it is defined 
% how each species contributes to Fred.
massBalanceElements.C = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0];
massBalanceElements.N = [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, cy1/cx1, cy2/cx2, 0, 0, cy1/cx1, 0, 0, 0, 0, 2, 0, 0];
massBalanceElements.Mn = [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
massBalanceElements.Fe = [0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0];
massBalanceElements.P = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, cz1/cx1, cz2/cx2, 1, 0, cz1/cx1, 0, 2, 0, 0, 0, 0, 0];
massBalanceElements.S = [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 2, cs1/cx1, cs2/cx2, 0, 0, cs1/cx1, 0, 0, 0, 0, 0, 0, 0];
massBalanceElements.Fred = [0, 0, 0, 0, 0, 2, 0.5, 0.25, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0]; 
 

 % R{i}            - the cell array of reactions. i=1:numRxns where numRxns 
 %                   are the total number of reactions in the network. 
 %                   For an irreversible reaction      aA + bB ---> cC + dD 
 %                   the R{i} written is 
 %                                                    a*A + b*B - c*C - d*D
 %                   For a reversible reaction         
 %                                                    aA + bB <---> cC + dD 
 %                   the R{i} will be
 %                                                    a*A + b*B - c*C - d*D  
 %                   and the R{i+1} will be
 %                                                    c*C + d*D - a*A - b*B
 % r{i}            - the cell array of reaction rates associated with each
 %                   reaction. It is a string variable
 %
 % isSolidReaction(i)- defines whether the reaction has the unit of a solid
 %                   i.e. umol/gr which equals to 1 or a solute with a unit 
 %                   of umol/cm3 which will equal to 0
 %
 % rName{i}         - cell array containing strings of names of reactions,
 %                    which are used for plotting purposes only

 
 % om1 oxidation by o2
 R{1}             = cx1*om1 + (cx1+0.25*cz1)*o2 - (cy1-2*cz1)*hco_3 -(cx1-cy1+2*cz1)*co_2 - (cx1-cy1+1.5*cz1)*h2o - cz1*hpo4 - cy1*nh4 - cs1*s2;
 r{1}             = 'komO * max(0,om1) * max(0,o2)/(o2 + kO) / cx1 * o2 / ( o2 + om1*1E-7)';
 isSolidReaction(1) = 1;
 rName{1}         = 'om1+o2->';
 
 % om1 denitrification 
 R{2}             = cx1*om1 + 0.8*cx1*no3 - (0.8*cx1+cy1-2*cz1)*hco_3 -(0.2*cx1-cy1+2*cz1)*co_2 - (0.6*cx1-cy1+2*cz1)*h2o - 0.4*cx1*n2 - cz1*hpo4 - cy1*nh4 - cs1*s2;
 r{2}             = 'komNO * max(0,om1) * max(0,no3)/(no3 + kNO) * kinO/(kinO + max(0,o2)) / cx1 * no3 / (no3+om1*1E-7)';
%r{3}             = 'komNO * max(0,om1) * max(0,no3)/(no3 + kNO) * kinO/(kinO + max(0,o2)) / cx1 * no3 / (no3+om1*3E-6) * om1 / (om1+no3*10)';
 isSolidReaction(2) = 1;
 rName{2}         = 'om1+no3->';

 % om1 oxidation by 1st MnO2 pool (mo_1)
 R{3}       	  = cx1*om1 + 2*cx1*mo_1 + (3*cx1+cy1-2*cz1)*co_2 + (cx1+cy1-2*cz1)*h2o - (4*cx1+cy1-2*cz1)*hco_3 - 2*cx1*mn - cz1*hpo4 - cy1*nh4 - cs1*s2;
 r{3}             = 'komMO * max(0,om1) * max(0,mo_1)/(abs(mo_1) + kMO) * kinO/(kinO + max(0,o2)) * kinNO/(kinNO + max(0,no3)) * mo_1 / (mo_1 + om1*1E-6) / cx1';
 isSolidReaction(3) = 1;
 rName{3}         = 'om1+mo_1->';

 % om1 oxidation by 1st Fe(OH) pool  (foh_1)
 R{4}             = cx1*om1 + 4*cx1*foh_1 + (7*cx1+cy1-2*cz1)*co_2 - (3*cx1-cy1+2*cz1)*h2o - (8*cx1+cy1-2*cz1)*hco_3 - 4*cx1*fe - cz1*hpo4 - cy1*nh4 - cs1*s2;
 r{4}             = 'komFOH * max(0,om1) * max(0,foh_1)/(foh_1 + kFOH) * kinO/(kinO + max(0,o2)) * kinNO/(kinNO + max(0,no3)) * kinMO/(kinMO + max(0,mo_1))/ cx1  * foh_1 / (foh_1 + om1*1E-6)';
 isSolidReaction(4) = 1;
 rName{4}         = 'om1+foh_1->';

 % om1 oxidation by sulfate
 R{5}             = cx1*om1 + 0.5*cx1*so4 + (cy1-2*cz1)*co_2 + (cy1-2*cz1)*h2o - (cx1+cy1-2*cz1)*hco_3 - 0.5*cx1*s2 - cz1*hpo4 - cy1*nh4 - cs1*s2;
 r{5}             = 'komSO * max(0,om1) * max(0,so4)/(so4 + kSO) * kinO/(kinO + max(0,o2)) * kinNO/(kinNO + max(0,no3)) * kinMO/(kinMO + max(0,mo_1)) * kinFOH/(kinFOH + max(0,foh_1)) /cx1 * so4 / (so4 + om1*1E-8)';
 isSolidReaction(5) = 1;
 rName{5}         = 'om1+so4->';

 % methanogenesis by om1
 R{6}             = cx1*om1 + (cy1-2*cz1)*h2o - (0.5*cx1-cy1+2*cz1)*co_2 - (cy1-2*cz1)*hco_3 - 0.5*cx1*ch4 - cz1*hpo4 - cy1*nh4 - cs1*s2;
 r{6}             = 'komCH * max(0,om1) * kinO/(kinO + max(0,o2)) * kinNO/(kinNO + max(0,no3)) * kinMO/(kinMO + max(0,mo_1)) * kinFOH/(kinFOH + max(0,foh_1))* kinSO/(kinSO + max(0,so4)) / cx1';
 isSolidReaction(6) = 1;
 rName{6}         = 'om1->ch4';
 
 % nitrification
 R{7}             = nh4 + 2*o2 + 2*hco_3 - no3 - 2*co_2 - 3*h2o;
 r{7}             = 'knhox * max(0,nh4) * max(0,o2) * o2 / ( o2 + nh4*1E-3)';
 isSolidReaction(7) = 0;
 rName{7}         = 'nh4+o2->';
 
 % Mn oxidation by O2
 R{8}             = 2*mn + o2 + 2*hco_3 - 2*mo_2 - 2*co_2 - 2*h2o;
 r{8}             = 'kmox * max(0,mn) * max(0,o2) * o2 / ( o2 + mn*1E-5)';
 isSolidReaction(8) = 0;
 rName{8}         = 'mn+o2->';
 
 % Fe oxidation by O2
 R{9}             = fe + 0.25*o2 + 2*hco_3 +0.5*h2o - 0.5*foh_1 - 0.5*foh_2 - 2*co_2;
 r{9}             = 'kfox * max(0,fe) * max(0,o2) * o2 / ( o2 + fe*1E-5)';
 isSolidReaction(9) = 0;
 rName{9}         = 'Fe+o2->';
 
 % H2S oxidation by O2
 R{10}             = s2 + 2*o2 + 2*hco_3 - so4 - 2*co_2 - 2*h2o;
 r{10}             = 'ksox * max(0,s2) * max(0,o2) * o2 / ( o2 + s2*1E-3)';
 isSolidReaction(10) = 0;
 rName{10}         = 's2+o2->';

 % Methane oxidation by O2
 R{11}             = ch4 + 2*o2 - co_2 - 2*h2o;
 r{11}             = 'kchox * max(0,ch4) * max(0,o2) * o2 / ( o2 + ch4*1E-2 ) * ch4 / ( ch4 + o2*1E-5)';
 isSolidReaction(11) = 0;
 rName{11}         = 'ch4+o2->';

 % Ammonium oxidation by 1st MnO2 pool(mo_1)
 R{12}             = 2*nh4 + 3*mo_1 + 5*h - 3*mn - n2 - 6*h2o;
 r{12}             = 'knhmo * max(0,nh4) * max(0,mo_1) * mo_1 / (mo_1 + nh4*1E-2 )';
 isSolidReaction(12) = 1;
 rName{12}         = 'nh4+mo_1->';
 
 % Ammonium oxidation by 2nd MnO2 pool(mo_2)
 R{13}             = 2*nh4 + 3*mo_2 + 5*h - 3*mn - n2 - 6*h2o;
 r{13}             = 'knhmx * max(0,nh4) * max(0,mo_2)  * mo_2 / (mo_2 + nh4*1E-2 )';
 isSolidReaction(13) = 1;
 rName{13}         = 'nh4+mo_2->';

 % Fe oxidation by 1st MnO2 pool(mo_1) 
 R{14}             = 2*fe + mo_1 + 2*hco_3 + 2*h2o - 1*foh_1 - 1*foh_2 - mn - 2*co_2;
 r{14}             = 'kfmo * max(0,mo_1) * max(0,fe) * mo_1 / (mo_1 + fe*1E1 )';
 isSolidReaction(14) = 1;
 rName{14}         = 'fe+mo_1->';
 
 % Fe oxidation by 2nd MnO2 pool(mo_2)
 R{15}             = 2*fe + mo_2 + 2*hco_3 + 2*h2o - 1*foh_1 - 1*foh_2 - mn - 2*co_2;
 r{15}             = 'kfmx * max(0,mo_2) * max(0,fe) * mo_2 / (mo_2 + fe*1E-8 )';
 isSolidReaction(15) = 1;
 rName{15}         = 'fe+mo_2->';
 
 % H2S oxidation by 1st MnO2 pool(mo_1)
 R{16}             = s2 + 4*mo_1 + 6*co_2 + 2*h2o - 4*mn - so4 - 6*hco_3;
 r{16}             = 'ksmo * max(0,mo_1) * max(0,s2) * mo_1 / (mo_1 + s2*1E-2 )';
 isSolidReaction(16) = 1;
 rName{16}         = 's2+mo_1->';
 
 % H2S oxidation by 2nd MnO2 pool(mo_2)
 R{17}             = s2 + 4*mo_2 + 6*co_2 + 2*h2o - 4*mn - so4 - 6*hco_3;
 r{17}             = 'ksmx * max(0,mo_2) * max(0,s2)';
 isSolidReaction(17) = 1;
 rName{17}         = 's2+mo_2->';
 
 % Fe(OH) 1st pool reduction by H2S
 R{18}             = s2 + 8*foh_1 + 14*co_2 - 8*fe - so4 - 14*hco_3 - 6*h2o;
 r{18}             = 'ksfo * max(0,foh_1) * max(0,s2) * foh_1 / (foh_1 + s2*1E1 )';
 isSolidReaction(18) = 1;
 rName{18}         = 's2+foh_1->';
 
  % Fe(OH) 2nd pool reduction by H2S
 R{19}             = s2 + 8*foh_2 + 14*co_2 - 8*fe - so4 - 14*hco_3 - 6*h2o;
 r{19}             = 'ksfx * max(0,foh_2) * max(0,s2) * foh_2 / (foh_2 + s2*1E-2 )';
 isSolidReaction(19) = 1;
 rName{19}         = 's2+foh_2->';
 
  % anaerobic methane oxidation
 R{20}             = ch4 + so4 + co_2 - 2*hco_3 -s2;
 r{20}             = 'kchso * max(0,ch4) * max(0,so4) * ch4 / ( ch4 + so4*1E-3 )';
 isSolidReaction(20) = 0;
 rName{20}         = 'ch4+so4->';
 
     % Mn oxidation by nitrate
 R{21}             = mn + 0.4*no3 + 1.6*hco_3 - mo_2 - 0.2*n2 - 1.6*co_2 - 0.8*h2o;
 r{21}             = 'kmoN * max(0,mn) * max(0,no3) * no3 / ( no3 + mn*1E-4)';
 isSolidReaction(21) = 0;
 rName{21}         = 'mn+no3->';
 
  % Ammonium oxidation by 1st FeOOH pool(foh_1)
 R{22}             = nh4 + 3*foh_1 + 5*h - 3*fe - 0.5*n2 - 6*h2o;
 r{22}             = 'knhfo * max(0,nh4) * max(0,foh_1) * foh_1 / (foh_1 + nh4*1E-8 )';
 isSolidReaction(22) = 1;
 rName{22}         = 'nh4+foh_1->';

  % Ammonium oxidation by 2nd FeOOH pool(foh_2)
 R{23}             = nh4 + 3*foh_2 + 5*h - 3*fe - 0.5*n2 - 6*h2o;
 r{23}             = 'knhfx * max(0,nh4) * max(0,foh_2) * foh_2 / (foh_2 + nh4*1E-8 )';
 isSolidReaction(23) = 1;
 rName{23}         = 'nh4+foh_2->';
 
  % iron sulfide precipitation
 R{24}             = s2 + fe - 2*h - fs;
 r{24}             = 'kIronSulfidePre * (max(0,fe)*max(0,s2)/(kEqIronSulfide*ph)-1) * (max(0,fe)*max(0,s2)/(kEqIronSulfide*ph^1)>=1)';
 isSolidReaction(24) = 1;
 rName{24}         = 's2+fe->';

 % iron sulfide dissolution
 R{25}             = fs + 2*h - fe - s2;
 r{25}             = 'kIronSulfideDis * max(0,fs) * (1-max(0,fe)*max(0,s2)/(kEqIronSulfide*ph^1)) * (max(0,fe)*max(0,s2)/(kEqIronSulfide*ph^1)<1) * fs / (fs + 1E-20 )';
 isSolidReaction(25) = 1;
 rName{25}         = 'fs->';

 % pyrite precipitation
 R{26}             = fs + s2 - pyrite + h2;
 r{26}             = 'kpyrpre * max(0,fs) * max(0,s2) * fs / (fs + s2*1E-20 )';
 isSolidReaction(26) = 1;
 rName{26}         = 'fs+s2->';
 
  % vivianite precipitation
 R{27}             = 3*fe + 2*hpo4 - vivianite - 2*h - 8*h2o;
 r{27}             = 'kvivpre * (((fe)^3*(hpo4)^2/kEqViv)^0.2 - 1) * ((max(0,fe))^3*(max(0,hpo4))^2/kEqViv>=1)';
 isSolidReaction(27) = 1;
 rName{27}         = 'fe+hpo4->';
 
 % vivianite dissolution
 R{28}             = vivianite + 2*h + 8*h2o - 3*fe - 2*hpo4;
 r{28}             = 'kvivdis * max(0,vivianite) * (1 - ((max(0,fe))^3*(max(0,hpo4))^2/kEqViv)^0.2) * ((max(0,fe))^3*(max(0,hpo4))^2/kEqViv<1)';
 isSolidReaction(28) = 1;
 rName{28}         = 'viv->';
 
  % vivianite dissolution by H2S
 R{29}             = 3*s2 + vivianite - 3*fs - 2*hpo4 - 4*h;
 r{29}             = 'ksviv * max(0,s2) * max(0,vivianite)';
 isSolidReaction(29) = 1;
 rName{29}         = 's2+viv->';

 % Manganese Carbonate precipitation
 R{30}             = mn + 2 * hco_3 - mo_3 - co_2 - h2o;
 r{30}             = 'kMnCarbonatePre * (mn*hco_3*kEqHCO3CO3/(ph*kEqMnCarbonate)-1) * (max(0,mn)*max(0,hco_3)*kEqHCO3CO3/(ph*kEqMnCarbonate)>=1) * mn / (mn + 1E-20 )';
 isSolidReaction(30) = 1;
 rName{30}         = 'mn+hco_3->';
 
  