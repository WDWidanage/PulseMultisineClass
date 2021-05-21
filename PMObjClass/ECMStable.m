function [Vlw, J] = ECMStable(theta, Il, OCV, timeVec, order)
% Stable equivalent cirucit model, without series capacitor
%
% Parameters are arranged as:
%      theta = Ro, Rp1,..,Rpn, tau1,...,taun
%
% W. D. Widanage 04/07/2014 (Serene)
%
dataPts = length(Il);       % Number of data points

% Check if weighting function is provided. else initialise to ones

% Change parameter vector to a column vector
theta = theta(:);

% Extract parameters 
Ro = theta(1);
RpAll = theta(2:1+order);
tauAll = theta(2+order:1+2*order);



% Initialise
Ip = zeros(order,1); 
Vl = zeros(dataPts,1);

% Initialise for Jacobian
J = zeros(dataPts,length(theta));
dIpdtau = zeros(order,1);

Vl(1) = OCV - Ro*Il(1) - RpAll'*Ip;
J(1,:) = [-Il(1), -Ip', -dIpdtau'];

for ii = 2:dataPts
    Ts = timeVec(ii)-timeVec(ii-1);                                    % Sampling interval
    expTau = exp(-Ts./tauAll);
    for jj = 1:order      
        Ip(jj,1) = expTau(jj)*(Ip(jj,1) - Il(ii)) + Il(ii);
    end
    
    Vl(ii) = OCV - Ro*Il(ii) - RpAll'*Ip;    
    
    % Create Jacobian matrix
    for jj = 1:order
        dIpdtau(jj,1) = expTau(jj)*(dIpdtau(jj,1) + Ip(jj,1)*Ts/tauAll(jj)^2 - Il(ii)*Ts/tauAll(jj)^2);
    end
    RdIpdtau = RpAll.*dIpdtau; 
    J(ii,:) = [-Il(ii), -Ip', -RdIpdtau'];
    

end

% Weight the output with the standard deviation over periods
Vl = Vl(:); % Change to a column vector
Vlw = Vl;
