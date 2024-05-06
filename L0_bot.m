% Script to output L0 bottle results for the statistical analysis.

clear; clc; close all;
addpath("baroneRoutines\");
addpath("func\");
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [3 3 28 15]);
% tmpT = "";

% Possible test cases.
principleAnalysis = false;  % main analysis
seasonalAnalysisKs = true;   % seasonality of statistics: K-S
seasonalAnalysisAd = true;  % seasonality of statistics: A-D
analyseStartYear = false;  % analyse effect (if any) of varying start year 
                            % on distributions
analyseEndYear = false;      % effect of varying end year
crn131 = false;             % analyse 2001-2021 data (to mirror CTD results)
nightAnalysis = false;      % analyse night-time 2001-2021 (to mirror CTD 
                            % results)
logAxes = true;            % output p-values as log values

if logAxes == true
    lp = "log/";
else
    lp = "";
end


%% Seasonal Analysis: K-S
if seasonalAnalysisKs == true

    pVals = [];
    
    % WINTER
    tmpT = "-01";

    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,1);
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % SPRING
    tmpT = "-02";
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,2);
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % SUMMER
    tmpT = "-03";
    
    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,3);
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % AUTUMN
    tmpT = "-04";

    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes,4);
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ks(2,[1 3])];

    season = [1 2 3 4];
    ax = figure;
    semilogy(season,pVals([1 18 35 52],:));
    lgd = legend("5","25",Location="south");
    lgd.Title.String = "pressure (dbar)";
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex",HandleVisibility="off"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex",HandleVisibility="off");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex",HandleVisibility="off"); hold off; title("chl $a$ (K-S)",Interpreter="latex");
    set(gca,"XTick",1:1:4,"XTickLabel",["winter","spring","summer","autumn"]);
    xlabel("season"); ylabel("$p$-value",Interpreter="latex");
    exportgraphics(ax,"figures/L0/bot/synthSsnl/chla_ks.png");

end

%% Seasonal Analysis: A-D
if seasonalAnalysisAd == true

    pVals = [];

    % WINTER
    tmpT = "-ad-01";

    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,1);
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];


    % SPRING
    tmpT = "-ad-02";

    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,pB,X,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,2);
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];


    % SUMMER
    tmpT = "-ad-03";
    
    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,3);
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];


    % AUTUMN
    tmpT = "-ad-04";
    
    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,"ad",logAxes,4);
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clear tmp ax; pVals = [pVals; ad([1 3])];

    season = [1 2 3 4];
    ax = figure;
    semilogy(season,pVals([1 18 35 52],:));
    lgd = legend("5","25",Location="south");
    lgd.Title.String = "pressure (dbar)";
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex",HandleVisibility="off"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex",HandleVisibility="off");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex",HandleVisibility="off"); hold off; title("chl $a$ (A-D)",Interpreter="latex");
    set(gca,"XTick",1:1:4,"XTickLabel",["winter","spring","summer","autumn"]);
    xlabel("season"); ylabel("$p$-value",Interpreter="latex");
    exportgraphics(ax,"figures/L0/bot/synthSsnl/chla_ad.png");
end

%% Principal Analysis
if principleAnalysis == true    
    % K-S
    tmpT = "";

    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1988-2021)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"chla" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Monovinyl Chl $a$","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Divinyl Chl $a$","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $b$","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $c1 + c2 + c3$","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: $\alpha$-carotene","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: 19' Butanoyloxyfucoxanthin","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Zeaxanthin","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Particulate Carbon","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Particulate Nitrogen","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Low-Level Phosphorus","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Low-Level Nitrogen","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Nitrate + Nitrite","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Phosphate","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Dissolved Oxygen","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Primary Production (Light-12)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % macrozooplankton
    tmp = importdata('data/L0/macrozoo_94-22_200.txt');
    [ax,~,~,pIn] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Macrozooplankton (1994-2021)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"macrozoo" + tmpT + ".png");
    %clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    %%%

    % A-D
    tmpT = "-ad";
    
    % chla
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,pB,X] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % mvchla 
    tmp = importdata('data/L0/mvchla_94-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Monovinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "mvchla" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % dvchla
    tmp = importdata('data/L0/dvchla_94-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Divinyl Chl $a$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "dvchla" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % chlb
    tmp = importdata('data/L0/chlb_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $b$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlb" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % chlc123
    tmp = importdata('data/L0/chlc123_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $c1 + c2 + c3$"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chlc123" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % alpha-caro
    tmp = importdata('data/L0/acar_94-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: $\alpha$-carotene"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "acar" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % but-19
    tmp = importdata('data/L0/but19_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: 19' Butanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "but19" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % hex-19
    tmp = importdata('data/L0/hex19_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: 19' Hexanoyloxyfucoxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "hex19" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % zeax
    tmp = importdata('data/L0/zeax_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Zeaxanthin"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "zeax" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % pc
    tmp = importdata('data/L0/pc_89-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Particulate Carbon"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pc" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % pn
    tmp = importdata('data/L0/pn_89-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Particulate Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "pn" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % llp
    tmp = importdata('data/L0/llp_88-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Low-Level Phosphorus"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "llp" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % lln
    tmp = importdata('data/L0/lln_88-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Low-Level Nitrogen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "lln" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % nit + nit
    tmp = importdata('data/L0/nitNit_88-21_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Nitrate + Nitrite"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "nit2" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % phos
    tmp = importdata('data/L0/phos_88-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Phosphate"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "phos" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % boxy
    tmp = importdata('data/L0/boxy_88-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Dissolved Oxygen"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "boxy" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % L-12 PP
    tmp = importdata('data/L0/l12_89-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Primary Production (Light-12)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "l12" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;

    % macrozooplankton
    tmp = importdata('data/L0/macrozoo_94-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ad",logAxes);
    sgtitle("L0: Macrozooplankton (1994-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/"+lp+"macrozoo" + tmpT + ".png");
    %clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
end

%% 2001-2021 chla (starting CRN 131, to match newer fluorometer)
if crn131 == true
    
    tmpT = "";
    
    % K-S
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2001-2021)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;
    
    % A-D
    tmpT = "-ad";
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    [ax,~,~] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2001-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;

else
    disp("Not analysing CRN 131 only data...");
end

%% Chl-a 2001-2022 NIGHT-TIME
if nightAnalysis == true
        
    tmpT = "";
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    
    hms = char(string(tmp.data(:,3)));
    mdy = char(string(tmp.data(:,2)));
    
    % test = "0" + t2(1,1:end-1);
    for i = 1:length(tmp.data)
        if hms(i,end) == " "
            hms(i,:) = "0" + hms(i,1:end-1);
        end
        if mdy(i,end) == " "
            mdy(i,:) = "0" + mdy(i,1:end-1);
        end
    end
    
    Y = double("20" + mdy(:,5:6));
    M = double("" + mdy(:,1:2));
    D = double("" + mdy(:,3:4));
    h = double("" + hms(:,1:2));
    m = double("" + hms(:,3:4));
    s = double("" + hms(:,5:6));
    
    T = datetime(Y,M,D,h,m,s);
    T2 = datetime(Y,M,D);
    T3 = datenum(T);
    
    Lat = 22.75;
    Lon = -158;
    [SunRiseSet,~,~,~,~,~] = suncycle(Lat,Lon,T2);
    for i = 1:length(SunRiseSet)
        tmp = SunRiseSet(i,:) - 10;
        for j = 1:2
            if tmp(j) < 0
                tmp(j) = tmp(j) + 24;
            end
        end
        rs(i,:) = tmp;
    end
    rs2 = hours(rs);
    rs2.Format = 'hh:mm';
    
    sunriseTime = hours(rs2(:,1)');
    sunsetTime = hours(rs2(:,2)');
    
    dayFrac = rem(T3,1);
    castTime = dayFrac*24;
    
    castAtNight = nan(length(T),1);
    
    for i = 1:length(T)
        if castTime(i) < sunriseTime(i) || castTime(i) > sunsetTime(i)
            castAtNight(i) = 1;
        end
    end
    
    nightCastIDs = [];
    
    for i=1:length(T)
        if(~isnan(castAtNight(i)))
            disp(i);
            nightCastIDs = [nightCastIDs i];
        end
    end
    
    tmp = importdata('data/L0/hplcChla_01-22_200.txt').data(nightCastIDs,:);
    [ax,~,~] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2001-2021, NIGHT)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21_night" + tmpT + ".png");
    clearvars -except tmpT analyseStartYear crn131 nightAnalysis logAxes lp;

else
    disp("Not analysing night-time only data...");
end

%% Year-by-year analysis for chla
% Here we move the start date of the analysis forward in time to see if
% the distribution of data has some dependence on time. It will only be
% checked if "analyseStartYear" is true.

if analyseStartYear==true
    
    % save p-values per year
    yearList = 1988:1:2016;
    pVals = [];

    % K-S
    tmpT = "";
    
    % 2016-2022
    tmp = importdata('data/L0/hplcChla_16-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2016-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_16-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2015-2022
    tmp = importdata('data/L0/hplcChla_15-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2015-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_15-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];    

    % 2013-2022
    tmp = importdata('data/L0/hplcChla_13-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2013-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_13-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2012-2022
    tmp = importdata('data/L0/hplcChla_12-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2012-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_12-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2011-2022
    tmp = importdata('data/L0/hplcChla_11-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2011-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_11-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2011-2022
    tmp = importdata('data/L0/hplcChla_11-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2011-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_11-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2010-2022
    tmp = importdata('data/L0/hplcChla_10-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2010-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_10-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2009-2022
    tmp = importdata('data/L0/hplcChla_09-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2009-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_09-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2008-2022
    tmp = importdata('data/L0/hplcChla_08-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2008-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_08-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2007-2022
    tmp = importdata('data/L0/hplcChla_07-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2007-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_07-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2006-2022
    tmp = importdata('data/L0/hplcChla_06-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2006-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_06-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2005-2022
    tmp = importdata('data/L0/hplcChla_05-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2005-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_05-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2004-2022
    tmp = importdata('data/L0/hplcChla_04-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2004-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_04-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2003-2022
    tmp = importdata('data/L0/hplcChla_03-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2003-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_03-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2002-2022
    tmp = importdata('data/L0/hplcChla_02-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2002-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_02-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 2001-2022 done above.
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2001-2021)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % 2000-2022
    tmp = importdata('data/L0/hplcChla_00-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (2000-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_00-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1999-2022
    tmp = importdata('data/L0/hplcChla_99-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1999-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_99-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1998-2022
    tmp = importdata('data/L0/hplcChla_98-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1998-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_98-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1997-2022
    tmp = importdata('data/L0/hplcChla_97-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1997-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_97-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1996-2022
    tmp = importdata('data/L0/hplcChla_96-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1996-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_96-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1995-2022
    tmp = importdata('data/L0/hplcChla_95-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1995-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_95-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1994-2022
    tmp = importdata('data/L0/hplcChla_94-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1994-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_94-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1993-2022
    tmp = importdata('data/L0/hplcChla_93-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1993-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_93-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1992-2022
    tmp = importdata('data/L0/hplcChla_92-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1992-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_92-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1991-2022
    tmp = importdata('data/L0/hplcChla_91-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1991-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_91-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1990-2022
    tmp = importdata('data/L0/hplcChla_90-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1990-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_90-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    % 1989-2022
    tmp = importdata('data/L0/hplcChla_89-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,"ks",logAxes);
    sgtitle("L0: Chl $a$ (1989-2022)","Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_89-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % 1988-2022
    tmp = importdata('data/L0/hplcChla_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_88-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];
    
    %% chla (A-D): go back year by year.
    tmpT = "-ad";
    
    % 2016-2022
    tmp = importdata('data/L0/hplcChla_16-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2016-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_16-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2015-2022
    tmp = importdata('data/L0/hplcChla_15-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2015-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_15-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2014-2022
    tmp = importdata('data/L0/hplcChla_14-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2014-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_14-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2013-2022
    tmp = importdata('data/L0/hplcChla_13-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2013-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_13-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2012-2022
    tmp = importdata('data/L0/hplcChla_12-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2012-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_12-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2011-2022
    tmp = importdata('data/L0/hplcChla_11-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2011-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_11-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2010-2022
    tmp = importdata('data/L0/hplcChla_10-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2010-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_10-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2009-2022
    tmp = importdata('data/L0/hplcChla_09-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2009-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_09-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2008-2022
    tmp = importdata('data/L0/hplcChla_08-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2008-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_08-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2007-2022
    tmp = importdata('data/L0/hplcChla_07-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2007-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_07-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2006-2022
    tmp = importdata('data/L0/hplcChla_06-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2006-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_06-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2005-2022
    tmp = importdata('data/L0/hplcChla_05-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2005-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_05-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2004-2022
    tmp = importdata('data/L0/hplcChla_04-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2004-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_04-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2003-2022
    tmp = importdata('data/L0/hplcChla_03-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2003-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_03-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2002-2022
    tmp = importdata('data/L0/hplcChla_02-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2002-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_02-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 2001-2022 done above.
    tmp = importdata('data/L0/hplcChla_01-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2001-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_01-21" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2000-2022
    tmp = importdata('data/L0/hplcChla_00-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (2000-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_00-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1999-2022
    tmp = importdata('data/L0/hplcChla_99-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1999-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_99-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1998-2022
    tmp = importdata('data/L0/hplcChla_98-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1998-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_98-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1997-2022
    tmp = importdata('data/L0/hplcChla_97-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1997-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_97-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1996-2022
    tmp = importdata('data/L0/hplcChla_96-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1996-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_96-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1995-2022
    tmp = importdata('data/L0/hplcChla_95-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1995-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_95-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1994-2022
    tmp = importdata('data/L0/hplcChla_94-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1994-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_94-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1993-2022
    tmp = importdata('data/L0/hplcChla_93-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1993-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_93-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1992-2022
    tmp = importdata('data/L0/hplcChla_92-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1992-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_92-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1991-2022
    tmp = importdata('data/L0/hplcChla_91-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1991-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_91-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1990-2022
    tmp = importdata('data/L0/hplcChla_90-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1990-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_90-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];
    
    % 1989-2022
    tmp = importdata('data/L0/hplcChla_89-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1989-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_89-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 1988-2022
    tmp = importdata('data/L0/hplcChla_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/" + lp + "chla_89-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    %% start-year analysis synthesis plot
    ax = figure;
    subplot(1,2,1)
    semilogy(flip(yearList),pVals(1:29,:)); grid on;
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex"); hold off; title("K-S"); xlim([1988 2016]);
    legend("5 dbar","25 dbar",Location="northwest"); xlabel("starting year"); ylabel("$p$-value",Interpreter="latex");
    subplot(1,2,2)
    semilogy(flip(yearList),pVals(30:58,:)); grid on;
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex"); hold off; title("A-D"); xlim([1988 2016]);
    legend("5 dbar","25 dbar",Location="northwest"); xlabel("starting year"); ylabel("$p$-value",Interpreter="latex");
    sgtitle("does lognormality of surface chl $a$ depend on the start year of the analysis?","Interpreter","latex")
    exportgraphics(ax,"figures/L0/bot/startYear/surfaceChla.png");


else
    disp("Not analysing year by year...");
end

%% Analyse by end year
if analyseEndYear

    yearList = 1994:1:2022;
    pVals = [];
    % K-S
    tmpT = "-ks";

    % -2022
    tmp = importdata('data/L0/hplcChla_88-22_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2021
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-21" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2020
    tmp = importdata('data/L0/hplcChla_88-20_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2020)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-20" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2019
    tmp = importdata('data/L0/hplcChla_88-19_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2019)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-19" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2018
    tmp = importdata('data/L0/hplcChla_88-18_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2018)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-18" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2017
    tmp = importdata('data/L0/hplcChla_88-17_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2017)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-17" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2016
    tmp = importdata('data/L0/hplcChla_88-16_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2016)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-16" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2015
    tmp = importdata('data/L0/hplcChla_88-15_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2015)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-15" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2014
    tmp = importdata('data/L0/hplcChla_88-14_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2014)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-14" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2013
    tmp = importdata('data/L0/hplcChla_88-13_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2013)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-13" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2012
    tmp = importdata('data/L0/hplcChla_88-12_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2012)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-12" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2011
    tmp = importdata('data/L0/hplcChla_88-11_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2011)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-11" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2010
    tmp = importdata('data/L0/hplcChla_88-10_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2010)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-10" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2009
    tmp = importdata('data/L0/hplcChla_88-09_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2009)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-09" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2008
    tmp = importdata('data/L0/hplcChla_88-08_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2008)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-08" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2007
    tmp = importdata('data/L0/hplcChla_88-07_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2007)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-07" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2006
    tmp = importdata('data/L0/hplcChla_88-06_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2006)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-06" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2005
    tmp = importdata('data/L0/hplcChla_88-05_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2005)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-05" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2004
    tmp = importdata('data/L0/hplcChla_88-04_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2004)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-04" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2003
    tmp = importdata('data/L0/hplcChla_88-03_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2003)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-03" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2002
    tmp = importdata('data/L0/hplcChla_88-02_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2002)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-02" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2001
    tmp = importdata('data/L0/hplcChla_88-01_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2001)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-01" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -2000
    tmp = importdata('data/L0/hplcChla_88-00_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-2000)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-00" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -1999
    tmp = importdata('data/L0/hplcChla_88-99_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-1999)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-99" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -1998
    tmp = importdata('data/L0/hplcChla_88-98_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-1998)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-98" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -1997
    tmp = importdata('data/L0/hplcChla_88-97_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-1997)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-97" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -1996
    tmp = importdata('data/L0/hplcChla_88-96_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-1996)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-96" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -1995
    tmp = importdata('data/L0/hplcChla_88-95_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-1995)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-95" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];

    % -1994
    tmp = importdata('data/L0/hplcChla_88-94_200.txt');
    [ax,ks] = L0_helper(tmp,50,'ks');
    sgtitle("L0: Chl $a$ (1988-1994)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_88-94" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ks(2,[1 3])];


    % A-D
    tmpT = "-ad";

    % -2022
    tmp = importdata('data/L0/hplcChla_88-22_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2022)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-22" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % -2021
    tmp = importdata('data/L0/hplcChla_88-21_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2021)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-21" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % -2020
    tmp = importdata('data/L0/hplcChla_88-20_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2020)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-20" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2019
    tmp = importdata('data/L0/hplcChla_88-19_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2019)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-19" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2018
    tmp = importdata('data/L0/hplcChla_88-18_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2018)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-18" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2017
    tmp = importdata('data/L0/hplcChla_88-17_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2017)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-17" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2016
    tmp = importdata('data/L0/hplcChla_88-16_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2016)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-16" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2015
    tmp = importdata('data/L0/hplcChla_88-15_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2015)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-15" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2014
    tmp = importdata('data/L0/hplcChla_88-14_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2014)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-14" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2013
    tmp = importdata('data/L0/hplcChla_88-13_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2013)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-13" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2012
    tmp = importdata('data/L0/hplcChla_88-12_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2012)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-12" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2011
    tmp = importdata('data/L0/hplcChla_88-11_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2011)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-11" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2010  
    tmp = importdata('data/L0/hplcChla_88-10_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2010)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-10" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2009  
    tmp = importdata('data/L0/hplcChla_88-09_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2009)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-09" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2008  
    tmp = importdata('data/L0/hplcChla_88-08_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2008)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-08" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2007  
    tmp = importdata('data/L0/hplcChla_88-07_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2007)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-07" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2006  
    tmp = importdata('data/L0/hplcChla_88-06_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2006)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-06" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2005  
    tmp = importdata('data/L0/hplcChla_88-05_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2005)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-05" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2004  
    tmp = importdata('data/L0/hplcChla_88-04_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2004)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-04" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2003  
    tmp = importdata('data/L0/hplcChla_88-03_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2003)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-03" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2002 
    tmp = importdata('data/L0/hplcChla_88-02_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2002)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-02" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2001 
    tmp = importdata('data/L0/hplcChla_88-01_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2001)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-01" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 2000
    tmp = importdata('data/L0/hplcChla_88-00_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-2000)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-00" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 1999
    tmp = importdata('data/L0/hplcChla_88-99_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-1999)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-99" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 1998
    tmp = importdata('data/L0/hplcChla_88-98_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-1998)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-98" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 1997
    tmp = importdata('data/L0/hplcChla_88-97_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-1997)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-97" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 1996
    tmp = importdata('data/L0/hplcChla_88-96_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-1996)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-96" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 1995
    tmp = importdata('data/L0/hplcChla_88-95_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-1995)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-95" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % 1994
    tmp = importdata('data/L0/hplcChla_88-94_200.txt');
    [ax,~,~,~,~,ad] = L0_helper(tmp,50,'ad');
    sgtitle("L0: Chl $a$ (1988-1994)"+tmpT,"Interpreter","latex");
    exportgraphics(ax,"figures/L0/bot/endYear/" + lp + "chla_89-94" + tmpT + ".png");
    clear ax tmp; pVals = [pVals; ad([1 3])];

    % end-year analysis synthesis plot
    ax = figure;
    subplot(1,2,1)
    semilogy(flip(yearList),pVals(1:length(pVals(:,1))/2,:)); grid on;
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex"); hold off; title("K-S"); xlim([yearList(1) yearList(end)]);
    legend("5 dbar","25 dbar",Location="northwest"); xlabel("ending year"); ylabel("$p$-value",Interpreter="latex");
    %set(gca,"XTick",1:1:length(yearList),"XTickLabel",yearList);

    subplot(1,2,2)
    semilogy(flip(yearList),pVals(length(pVals(:,1))/2 + 1:end,:)); grid on;
    hold on; yline(0.05,"-","$\alpha$ = 0.05",Interpreter="latex"); yline(0.005,'--',"$\alpha$ = 0.005",Interpreter="latex");
    yline(0.1,'--',"$\alpha$ = 0.1",Interpreter="latex"); hold off; title("A-D"); xlim([yearList(1) yearList(end)]);
    legend("5 dbar","25 dbar",Location="northwest"); xlabel("ending year"); ylabel("$p$-value",Interpreter="latex");
    sgtitle("does lognormality of surface chl $a$ depend on the end year of the analysis?","Interpreter","latex")
    %set(gca,"XTick",1:1:length(yearList),"XTickLabel",yearList);
    exportgraphics(ax,"figures/L0/bot/endYear/"+lp+"surfaceChla.png");
end