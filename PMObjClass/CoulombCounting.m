function [soc, ampSec] = CoulombCounting(I,time,Cn_c,Cn_d,SoC0)
% Coulomb counting to estimate state-of-charge
% Perform current integration and normalise with resepect to battery
% capacity. Cn_c and Cn_d are the measured battery capacity when charging
% and discharging respectively at a given temperature.
%
% Integration is performed using trapezoidal method with saturation limits
% of 0 and 1
%
% Input arguments:
%   I: Current vector (A), size N x 1
%   time: time vector (s), size N x 1
%   Cn_c: Cell capacity when charging (Ah), size 1 x 1
%   Cn_d: Cell capacity when discharging (Ah), size 1 x 1
%   SoC0: Initial state-of-charge (%), size 1 x 1
%
% Output arguments:
% soc: Remaining state-of-charge (%), size N x 1
% ampSec: Ampere seconds (As), size N x 1
%
% W.D. Widanage 09/02/2013 (Boogie woogie!)

N = length(I);

% Initialise
soc = zeros(N,1);
ampSec = zeros(N,1);
soc(1) = SoC0;
ampSec(1) = 0;
inc = 0;

for ii = 1:N-1
    
    delT = time(ii+1)-time(ii);
    if I(ii)>=0
        inc = ((I(ii)+I(ii+1))*delT/2)/(Cn_d*3600);  % trapezoid increment approximation method and normalise with capacity (As)
    elseif I(ii)<0
        inc = ((I(ii)+I(ii+1))*delT/2)/(Cn_c*3600);  % trapezoid increment approximation method and normalise with capacity (As)
    end
    
    % Perform integration with saturation limits of 1 and 0
    if  soc(ii)-inc >1
        soc(ii+1) = 1;
    elseif soc(ii)-inc <0
        soc(ii+1) = 0;
    else
        soc(ii+1) = soc(ii)-inc;    % Accumulate remaining capacity
    end
    ampSec(ii+1) = ampSec(ii) + (I(ii)+I(ii+1))*delT/2;
end
