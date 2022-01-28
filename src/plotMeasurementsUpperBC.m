% Plots measurements that constrain the upperBC. This is not part of the
% MEDIALAB package, but only serves for visualization purposes of data
% stored in Solid_flux.mat and Solute_bc.mat

clear all
load('Solid_flux.mat')

% Plot solid fluxes
figure
set(gcf,'Position',get(0,'ScreenSize'));                          %maximize the figure
subplot(5,2,1)
plot(Ba13_Bulk(:,1),Ba13_Bulk(:,2),'k')
title('Bulk 2013', 'FontWeight','bold')
subplot(5,2,3)
plot(Ba13_TOC(:,1),Ba13_TOC(:,2),'b')
title('TOC 2013', 'FontWeight','bold')
subplot(5,2,5)
plot(Ba14_TOC(:,1),Ba14_TOC(:,2),'--b')
title('TOC 2014', 'FontWeight','bold')
ylabel('Absolute fluxes', 'FontWeight','bold')
subplot(5,2,7)
plot(Ba13_Fe(:,1),Ba13_Fe(:,2),'r')
title('Fe 2013', 'FontWeight','bold')
subplot(5,2,9)
plot(Ba13_Mn(:,1),Ba13_Mn(:,2),'g')
title('Mn 2013', 'FontWeight','bold')
xlabel('Time', 'FontWeight','bold')

% Plot solid fluxes relative to Bulk flux
subplot(5,2,4)
plot(Ba13_TOC(:,1),Ba13_TOC(:,2)./interp1(Ba13_Bulk(:,1),Ba13_Bulk(:,2),Ba13_TOC(:,1)),'b')
title('TOC 2013', 'FontWeight','bold')
subplot(5,2,6)
plot(Ba14_TOC(:,1),Ba14_TOC(:,2)./interp1(Ba13_Bulk(:,1),Ba13_Bulk(:,2),Ba14_TOC(:,1)),'--b')
title('TOC 2014', 'FontWeight','bold')
ylabel('Fluxes Relative to Bulk Flux', 'FontWeight','bold')
subplot(5,2,8)
plot(Ba13_Fe(:,1),Ba13_Fe(:,2)./interp1(Ba13_Bulk(:,1),Ba13_Bulk(:,2),Ba13_Fe(:,1)),'r')
title('FeT 2013', 'FontWeight','bold')
subplot(5,2,10)
plot(Ba13_Mn(:,1),Ba13_Mn(:,2)./interp1(Ba13_Bulk(:,1),Ba13_Bulk(:,2),Ba13_Mn(:,1)),'g')
title('MnT 2013', 'FontWeight','bold')
xlabel('Time', 'FontWeight','bold')


load('Solute_bc.mat')
% Plot solute BC
figure
set(gcf,'Position',get(0,'ScreenSize').*[1 1 0.5 1]);                          %maximize the figure
subplot(6,1,1)
plot(Ba01_febc(:,1),Ba01_febc(:,2),'.k')
title('Fe 2013', 'FontWeight','bold')
subplot(6,1,2)
plot(Ba01_hp04bc(:,1),Ba01_hp04bc(:,2),'.b')
title('HPO_4 2013', 'FontWeight','bold')
subplot(6,1,3)
ylabel('Solute Concentrations', 'FontWeight','bold')
plot(Ba01_mnbc(:,1),Ba01_mnbc(:,2),'.r')
title('Mn 2013', 'FontWeight','bold')
subplot(6,1,4)
plot(Ba01_nh4bc(:,1),Ba01_nh4bc(:,2),'.g')
title('NH_4 2013', 'FontWeight','bold')
subplot(6,1,5)
plot(Ba01_no3bc(:,1),Ba01_no3bc(:,2),'.c')
title('NO_3 2013', 'FontWeight','bold')
subplot(6,1,6)
plot(Ba01_so4bc(:,1),Ba01_so4bc(:,2),'.m')
title('SO4 2013', 'FontWeight','bold')
xlabel('Time', 'FontWeight','bold')



