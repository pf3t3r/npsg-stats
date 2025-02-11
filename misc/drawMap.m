% Draw map of Station ALOHA.

clc; clear; close all;

land = readgeotable("landareas.shp");
rivers = readgeotable("worldrivers.shp");
cities = readgeotable("worldcities.shp");
load coastlines
latlim = [17.5 25]; lonlim = [-162.5 -152.5];

ax = worldmap(latlim,lonlim);
plotm(coastlat,coastlon);
hold on
bordersm
geoshow(ax,land,"FaceColor",[0.5 0.7 0.5])
h(1) = geoshow(22.75, -158, 'DisplayType', 'Point', 'Marker', 'o','MarkerSize',12,MarkerFaceColor="#1f78b4",MarkerEdgeColor="k");
h(2) = geoshow(ax,cities,MarkerSize=10,MarkerEdgeColor='k',MarkerFaceColor=[0.8 0.8 0.8],Marker="square");
hold off
legend([h(1),h(2)],'Station ALOHA','Honolulu',fontsize=15);
% title("Hawai'i");

%% NL: NWO-I
latlim2 = [50.7 53.7];
lonlim2 = [3 7.4];

figure;
ax = worldmap(latlim2,lonlim2);
plotm(coastlat,coastlon,HandleVisibility="off");
hold on
geoshow(ax,land,"FaceColor",[1 1 1],HandleVisibility="off")
bordersm
h(1) = geoshow(52.35673045941863, 4.949037037824388,'DisplayName','AMOLF - ARCNL - CWI - NIKHEF','DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#a6cee3",MarkerEdgeColor="k");
h(2) = geoshow(52.81350697499335, 6.396063112749066,'DisplayName','ASTRON','DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#1f78b4",MarkerEdgeColor="k");
h(3) = geoshow(51.44883174206774, 5.494747224743309, 'DisplayName','DIFFER','DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#b2df8a",MarkerEdgeColor="k");
h(4) = geoshow(53.00305725888978, 4.786915772094974,'DisplayName','NIOZ Texel', 'DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#33a02c",MarkerEdgeColor="k");
h(5) = geoshow(51.48874010670033, 4.056787780243237, 'DisplayName','NIOZ Yerseke','DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#fb9a99",MarkerEdgeColor="k");
h(6) = geoshow(52.333629170832474, 4.868229830263699, 'DisplayName','NSCR','DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#e31a1c",MarkerEdgeColor="k");
h(7) = geoshow(53.2406132852035, 6.534953491813899, 'DisplayName','SRON Groningen','DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#fdbf6f",MarkerEdgeColor="k");
h(8) = geoshow(52.16957434458291, 4.457959663765497, 'DisplayName','SRON Leiden','DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#ff7f00",MarkerEdgeColor="k");
% h(9) = geoshow(52.35706663422197, 4.94783365090008, 'DisplayName','ARCNL','DisplayType', 'Point', 'Marker', '.','MarkerSize',6,MarkerFaceColor="#fdbf6f",MarkerEdgeColor="k");
% h(10) = geoshow(52.356682138840895, 4.95197439692403, 'DisplayName','CWI','DisplayType', 'Point', 'Marker', '*','MarkerSize',6,MarkerFaceColor="#ff7f00",MarkerEdgeColor="k");
% h(11) = geoshow(52.3563414264984, 4.951234073496149, 'DisplayName','NIKHEF','DisplayType', 'Point', 'Marker', 'o','MarkerSize',6,MarkerFaceColor="#cab2d6",MarkerEdgeColor="k");
hold off
legend([h(1),h(2),h(3),h(4),h(5),h(6),h(7),h(8)],sprintf('AMOLF, ARCNL,\nCWI, NIKHEF'),'ASTRON','DIFFER','NIOZ Texel','NIOZ Yerseke','NSCR','SRON Groningen','SRON Leiden');