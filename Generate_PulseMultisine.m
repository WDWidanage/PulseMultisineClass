% Script to generate pulse multisines for battery cell characterisation and
% modelling.
%
% The first key information required when designing a pulse-multisine is the maximum
% applicable 10 s charge (Ccmax) and discharge (Cdmax) rates. This can be
% found in the data sheet. These values are used to ensure that the
% pulse-multisine signal does not exceed battery safe c-rate region.
%
% A pulse-multisine signal is a combination of a base-signal and a
% multisine (see figure below). 
%
%           Base-signal          +         Multisine
%
%  --------Ccmax
%     
%     |--|                                   fMax
%     |  |              T4                |-|
%     |  |---------|  |-----     +    ---/   \-\   /----
%      T1    T2    |--| aplha                   |-|
%
%  ---------Cdmax
%
% It has five design parameters which shapes the base-signal and frequency 
% bandwidth of the multisine.
%
%   T1: Time interval of largest base-signal pulse. A typical value
%      of 10 s is sufficient when Cdmax = Ccmax. If Cdmax ~= Ccmax a value
%      of 5 s can be used as a rule of thumb.
%   T2: Time interval of first base-signal rest period A typical value 
%       of 20 s is sufficient for the design. 
%   T4 :Time interval of second base-signal rest period. A typical value 
%       of 20 s is sufficient for the design. 
%   alpha: Fraction of the smaller base-signal pulse compared to maximum
%          allowed c-rate. A typical value of 0.6 is sufficient for the 
%          design. If Cdmax ~= Ccmax a value of 0.5 can be used as a rule 
%          of thumb.
%   fMax : Highest excited frequency in the multisine signal. This should 
%          be set to cover the frequency of interest. 1Hz is sufficient for
%          a drive-cycle.
%
% Setting the design parameters results in a signal with characteristics
% similar to a drive cycle in terms of the exicted power spectrum. The
% subsequent equivalent circuit model estimated will therefore be 
% parameterised more accuratley within the domain of battery use. For more
% information see reference below.
%
% The script generates one period of a pulse-multisine signal. Once the 
% signal is generated apply P periods of the signal to the battery and 
% measure the P voltage responses. A value of P = 5 is sufficient to 
% characterise and model the battery dyanmics.
% 
% Reference:
%  Widanage, W. D., Barai, A., Chouchelamane, G.H., Uddin, K., McGordon, 
%  A., Marco, J. and Jennings, P., "Design and use of multisine signals for 
%  Li-ion battery equivalent circuit modelling. Part 1: Signal design", 
%  Journal of Power Sourcers, 324, pp. 70-78. 
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, p.refSig.K. 28/06/2016 
% W.D. Widanage 28/06/2016 (The Call of Ktulu)

clear
close all
addpath('PMObjClass')
%% Design parameters. Only this block of code needs to change for different batteries and different maximum discharge and charge rates

refCell = 'NCA_3Ah';
cRate = 3.03;       % C-rate of the battery
cDmax = 3;          % Maximum applicable 10 s discharge c-rate [-], for a given SoC and temperature, as specified by the manufacture
cCmax = 3;          % Maximum applicable 10 s dharge c-rate [-], for a given SoC and temperature, as specified by the manufacture

T1 = 10;            % Time interval [s] of the largest base-signal pulse
T2 = 20;            % Time interval [s] of the first base-signal rest period
T4 = 20;            % Time interval [s] of second base-signal rest period
alpha = 0.6;        % Fraction [-] of smaller base-signal pulse compared to maximum allowed c-rate
fMax = 1;           % Highest excited frequency [Hz] in the multisine signal. This should be set to cover the frequency of interest. 1Hz is sufficient for a drive-cycle
fs = 10;            % Sampling frequency [Hz] at which the cell cylcer is run. Normmaly it is at 10Hz. Note that fs should always be fs >= 2*fMax.

% Used for filename saving. No influence on signal design
refSoC = 50;           % SoC [-] at which pulse-multisine is expected to be applied
refTemp = 25;          % Temperature [degC] at which pulse-multisine is expected to be applied
saveSignal = 'n';      % If saveSignal is 'y', the script will: 
                       %    - Generate a text file to import the signal to a cell cycler (the format of text file will need adapting for the cycler of interest). 
                       %      This saved text file signal is one period of the pulse-multisien with time and current as the two columns. 
                       %    - Save a PMObj variable 'p' which you will need to load once the 
                       %      pulse-multisine signal is applied and the voltage
                       %      response is measured. See "estimateNLECM.m" on how
                       %      this PMObj is loaded and the measured current and
                       %      voltage are used to then estimate a NL-ECM. 

%% Create pulse-multisine signal
sigProp = struct('cRate', cRate, 'cDmax',cDmax,'cCmax',cCmax,'T1',T1,'T2',T2,'T4',T4,'alpha',alpha,'fMax',fMax,'fs',fs,'refCell',refCell,'refSoC',refSoC,'refTemp',refTemp);       % Create a structure variable with the signal parameters
p = PMObj(sigProp);


%% Plots
maxIrate = max(p.refSig.pmSignal)/cRate;
minIrate = min(p.refSig.pmSignal)/cRate;
currRms = rms(p.refSig.pmSignal);
BS = fft(p.refSig.baseSignal);

freqRec = [0:p.refSig.N-1]/(p.refSig.N)*fs;

U = fft(p.refSig.pmSignal);
excitedLines = (p.refSig.harmExc)+1;        % Excited lines realtive to record
freqExc = (excitedLines-1)*fs/(p.refSig.N);
freqSupp = p.refSig.harmSupp*fs/p.refSig.N;


figure()
plot(p.refSig.timeVec,p.refSig.pmSignal,'- .');
xlabel('Time (s)'); ylabel ('Current (A)'); title(['Pulse-multisine. Cmin: ',num2str(minIrate),' Cmax: ',num2str(maxIrate), ' Rms: ',num2str(currRms)])

figure()
hist(p.refSig.pmSignal,18)
xlabel('Signal amplitude'); ylabel('Frequency of occurrence'); title('Amplitude distribution of pulse-multisine.')

% Perform Coulomb counting to check SoC variation
[soc, ampSec] = CoulombCounting(p.refSig.pmSignal,p.refSig.timeVec,cRate,cRate,refSoC/100);
socChange = [max(soc)-min(soc)]*100;

figure()
plot(p.refSig.timeVec,soc*100)
xlabel('Time (s)'); ylabel('SoC (%)'); title(['SoC Change: ', num2str(socChange),'%'])

figure();
plot(freqRec,abs(BS),'-x',freqExc,abs(U(excitedLines)),'-or',freqSupp,abs(U(p.refSig.harmSupp+1)),'go');
xlabel('Frequency (Hz)')
ylabel('FFT magnitude (abs)')
legend('Base-signal fft','Pulse-multisine fft','Suppressed harmonics')
title('Amplitude spectrum of pulse-multisine')
xlim([0, fMax ])

%% Save signal
if ismember(saveSignal,['Y','y'])
    % Reference signal with harmonic specification
    save([refCell,'_PMS_',num2str(refSoC),'per_',num2str(refTemp),'degC'],'p')
    
    
    % Save for bitrode as a text file
    textFileName = [refCell,'_RefCurr_',num2str(refSoC),'per_',num2str(refTemp),'degC.txt'];
    currBitrode = [0;-p.refSig.pmSignal];
    timeBitrode = [0,1:p.refSig.N]'/fs;
    dataBitrode = [timeBitrode,currBitrode];
    save(textFileName,'dataBitrode','-ascii','-tabs')
end