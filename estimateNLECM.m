% Test PMObj class: This script demostrates on how to estimate a non-linear
% ECM (NLECM) based on a generated pulse-multisine signal using  the PMObj
% class.
%
% This requires that the reference pulse-multisine is saved. This is tehn
% loaded and the measured data is passed into the reference pulse-multisien
% object to estimate the NLECM.
%
% Reference:
%  Widanage, W. D., Barai, A., Chouchelamane, G.H., Uddin, K., McGordon, 
%  A., Marco, J. and Jennings, P., "Design and use of multisine signals for 
%  Li-ion battery equivalent circuit modelling. Part 2: Signal design", 
%  Journal of Power Sourcers, 324, pp. 70-78. 
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 01/08/2019 (Suffering vs satisfaction)

clear
clear class
clc
close all
addpath('PMObjClass')

%% Esimate the NL-ECM model parameters

% Once measurements are avaialble load the current and voltage response
% back into the pulse-multisine object that was created with the 
% "Generate_PulseMultisine.m" script. 


load('Reference_PMS_50per_25degC.mat') % Load the reference pulse-multisine signal that was generated with the "Generate_PulseMultisine.m" script. This is a PMObj type
load('Measured_PMS_50per_25degC.mat')  % Load the measured data

p.measSig.measCurr = -data.Current;    % Assign the measured pulse multisine current, -ve is charge
p.measSig.measVol = data.Voltage;      % Assign the meaasured voltage response
p.measSig.measTime = data.TotalTime;   % Assign the measured time
p.measSig.P = 5;

p.estSetting.ECMOrder = 1;             % State the ECM order to esitmate
p.estSetting.fCutOff = 0.7;            % Set any cut-off frequency

p.estNLECM;                            % Estimate the NLECM

p.plotAll;                             % Plot all the figures if required for analysis

%% Print the NL-ECM parameters for convenience

% Parameters related to the linear ECM model
Ro = p.nlECMPara.Ro     % The series Ohmic resistance [Ohms]
Rp = p.nlECMPara.Rp     % The polarisation resistances [Ohms]
Tau = p.nlECMPara.Tau   % The time constants [s]

% Parameters related to the non-linear sigmoid function. See equation (11)
% in the reference mentioned above
c1 = p.nlECMPara.c1    % Sigmoid fucntion coefficient [-]
c2 = p.nlECMPara.c2    % Sigmoid function coefficient [-]
