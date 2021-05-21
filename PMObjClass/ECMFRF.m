function [ecmFRFSplit, JSplit] = ECMFRF(theta,w,varargin)
% Return the frequency response of a stable ECM (without a series
% capacitor). 
%
% Inputs:
%       w: Angular frequency vector (rad/sample) size N x 1
%   theta: Parameter vector, [R0, Rp1,...,RpN, tau1,...,tauN], size nTheta x 1
%
% Outputs:
%   ecmFRF: FRF of ECM, size N x 1
%        J: Jacobian N x nTheta
%
% Optional input arguments. Create a structure variable with the following fields:
%   fs: Sampling frequency (Hz). Defualt fs = 1Hz
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 06/06/2016 (Waking up!)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only


p = inputParser; % Create an input parse object to handle positional and property-value arguments

% Create variable names and assign default values after checking the value
addRequired(p,'theta', @isnumeric);
addRequired(p,'w', @isnumeric);


% Optional parameteres
addParameter(p,'fs',1)

% Re-parse parObj
parse(p,theta,w,varargin{:})

% Caluclate FRF
order = (length(p.Results.theta)-1)/2;
R0 = p.Results.theta(1);
Rp = p.Results.theta(2:order+1);
tau = p.Results.theta(order+2:end);

% Angular frequency
w = p.Results.w * p.Results.fs;

for oo = 1:order
    rcFRF(:,oo) = Rp(oo)./(tau(oo)*1i*w + 1);
    
    dZdRp(:,oo) = 1./(tau(oo)*1i*w +1);
    dZdTau(:,oo) = -Rp(oo)./((tau(oo)*1i*w +1).^2)*1i.*w;
end

ecmFRF = R0 + sum(rcFRF,2);

ecmFRFSplit = [real(ecmFRF);imag(ecmFRF)];

% Jacobian matrix
J = [ones(length(w),1), dZdRp, dZdTau];
JSplit = [real(J);imag(J)];