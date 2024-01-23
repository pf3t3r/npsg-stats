clc;clear;close all;
% 
land = readgeotable("landareas.shp");
rivers = readgeotable("worldrivers.shp");
cities = readgeotable("worldcities.shp");

%%

latlim = [17.5 25];
lonlim = [-162.5 -152.5];
ax = worldmap(latlim,lonlim);

% worldmap("Pacific",[10 30]);
load coastlines
plotm(coastlat,coastlon);
hold on
geoshow(ax,land,"FaceColor",[0.5 0.7 0.5])
h(1) = geoshow(22.75, -158, 'DisplayType', 'Point', 'Marker', 'o', 'Color', 'red');
h(2) = geoshow(ax,cities);
hold off
legend([h(1),h(2)],'Station ALOHA','Honolulu');
% title("Hawai'i");