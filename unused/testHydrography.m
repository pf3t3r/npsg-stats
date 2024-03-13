close all;clc;clear;

sdata = importdata("data\L0\sp_88-21_4800.txt").data;
tdata = importdata("data\L0\t_88-21_4800.txt").data;

p = sdata(:,3); sp = sdata(:,4); t = tdata(:,4);
sp(sp==-9) = nan; t(t==-9) = nan;
long = -158;
lat = 22.75;

sa = gsw_SA_from_SP(sp,p,long,lat);
ct = gsw_CT_from_t(sa,t,p);
z = gsw_z_from_p(p,lat);

%%
figure
scatter(sa,p,[],[0.6 0.6 0.6],'.');
set(gca,'YDir','reverse');
ylabel('P [dbar]'); xlabel('$S_A$ [g/kg]','Interpreter','latex');

%%
figure
scatter(ct,p,[],[0.6 0.6 0.6],'.');
set(gca,'YDir','reverse');
ylabel('P [dbar]'); xlabel('$\Theta$ [$^\circ$]','Interpreter','latex');

%%
gamma_n = eos80_legacy_gamma_n(sp,t,p,long,lat);
sigma0 = gsw_sigma0(sa,ct);
sigma4 = gsw_sigma4(sa,ct);

%%
figure
scatter(gamma_n,p,[],[0.4 0.4 0.4],'.'); hold on
scatter(sigma0,p,[],'k','.');
scatter(sigma4,p,[],'r','o'); hold off
set(gca,'YDir','reverse'); legend('$\gamma^N$','$\sigma_{\theta}$','$\sigma_4$','Interpreter','latex');
ylabel('P [dbar]'); xlabel('Density [kg m$^{-3}$]','Interpreter','latex');

%%
[Nsquared,pmid] = gsw_Nsquared(sa,ct,p,lat);

figure
scatter(Nsquared,pmid,[],[0.4 0.4 0.4],'.');
set(gca,'YDir','reverse'); xlim([-1e-3 1e-3]);
ylim([50 180]);
ylabel('P [dbar]'); xlabel('N$^2$ [rad$^2 s^{-2}$]','Interpreter','latex');

%%
% Omega = 7.292e-5;
% fC = 2*Omega*sin(deg2rad(lat));
% PV = nan(length(p),1);
% 
% DF = gradient(sigma0);
% PV = fC.*(sigma0 + 1000).*DF;
% 
% figure
% scatter(PV,p,[],[0.4 0.4 0.4],'.'); set(gca,'YDir','reverse');
% xlim([-1e-3 1e-3]);
% ylabel('P [dbar]'); xlabel('PV [m$^{-1}s^{-1}$]','Interpreter','latex');