clear;clc;close all;
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 5 10 10]);

data = importdata('data/hots-chloropigment.txt').data;

p = data(:,3);
f = data(:,4);

f(f==-9) =nan;

pp = reshape(p,[101 329]);
ff = reshape(f,[101 329]);
ffm = mean(ff,2,"omitnan");
f5 = prctile(ff',5);
f95 = prctile(ff',95);

ax = figure;
plot(f,p,'.',Color=[0.8 0.8 0.8],DisplayName="raw data");
hold on
plot(ffm,pp(:,1),'-',"Color",[0 0 0],DisplayName="mean");
plot(f5,pp(:,1),'-',"Color",[0.5 0.5 0.5],DisplayName="5%");
plot(f95,pp(:,1),'-',"Color",[0.5 0.5 0.5],DisplayName="95%");
hold off
set(gca,"YDir","reverse");
legend();
% title("L0 Fluorescence 1988-2022",Interpreter="latex");
xlabel("chl-$a$ fluorescence [$\mu$g/L]",Interpreter="latex");
% ylabel("Pressure [dbar]",Interpreter="latex");
% yticklabels({});
yticks(0:20:200);
exportgraphics(ax,'figures/L0/meanProfile.png');
