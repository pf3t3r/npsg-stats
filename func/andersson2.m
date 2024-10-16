% Equation 2 from Andersson (2021). This equation describes the growth of
% phytoplankton. It assumes that growth comprises and deterministic
% constant component 'muk' and a stochastic component 'sigk*eta', where
% sigk is the magnitude of the stochastic variation and eta is the
% time-varying rate. For our analysis we introduce an additional mixing
% term 'M'. This can be constant or also time-varying. This equation is
% designed to be used with the ode45() solver.

function dxdt = andersson2(x,muk,sigk,M,timeVary)
if nargin < 5
    timeVary = false;
end
if nargin<4
    eta = rand();
    dxdt = -(muk + sigk*eta)*x;
else
    eta = rand();
    if timeVary == true
        dxdt = -(muk + sigk*eta)*x + M*rand();
    else
        dxdt = -(muk + sigk*eta)*x + M;
    end
end
end

