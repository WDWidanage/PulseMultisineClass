function [u] = PulseMultisine(cDmax, cCmax, cRate, fs, T1, T2, T4, alpha, fMax)
%
% Design a pulse-multisine signal for a given a maximium applicable
% discharge and charge current pulse.
%
% Mandotory input argumetns
%   cDmax: Maximum applicable discharge pulse, a double size 1 x 1 (usually obtianed from data sheet)
%   cCmax: Maximum applicable charge pulse, a double size 1 x 1 (usually obtained from data sheet)
%   cRate: The C-rate of the battery
%
% Optinal argumnets:
%   fs: Sampling frequency, default fs = 1Hz, double size 1 x 1
%   T1: Pulse length of largest applicable pulse, default T1 = 10s, double size 1 x 1
%   T2: First rest preiod of base-signal, default T2 = 20s, double size 1 x 1
%   T4: Second rest period of base-signal, default T4 = 20s, double size 1 x 1
%   alpha: Sclae factor of largest pulse, default alpha = 0.6, a double 0 < alpha < 1, size 1 x 1
%   fMax: Maximum frequency of multisine signal, default fMax = 1, double size 1 x 1
%
% Output arguments:
%   u: Pulse-multisine and time vector and pulse-multisine properties for reference, a 17 x 1 structure variable
%
% Reference:
%  Widanage, W. D., Barai, A., Chouchelamane, G.H., Uddin, K., McGordon, A.,
%  Marco, J. and Jennings, P., "Design and use of multisine signals for Li-ion battery equivalent circuit modelling. Part 1: Signal design", 
%  Journal of Power Sourcers, 324, pp. 70-78. 
%
% Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 26/04/2015 (Hope)
% All Rights Reserved
% Software may be used freely for non-comercial purposes only


if nargin < 4 || isempty(fs)
    fs = 1;
end
if nargin < 5 || isempty(T1)
    T1 = 10;
end
if nargin < 6 || isempty(T2)
    T2 = 20;
end
if nargin < 7 ||isempty(T4)
    T4 = 20;
end
if nargin < 8 || isempty(alpha)
    alpha = 0.6;
end
if nargin < 9 || isempty(fMax)
    fMax = 1;
end


% Cmin and Cmax
cMin = min([cDmax,cCmax]);
cMax = max([cDmax,cCmax]);

% C-rate of smallest pulse
C2 = alpha*cMin;

% beta, scale factor of multisine signal
beta = 1 - alpha;

% gamma, scale factor of largest pulse
gamma = (cMax - beta*cMin)/cMax;

% C-rate of largest pulse
C1 = gamma * cMax;

% T3, time length of smallest pulse
T3 = C1*T1/C2;

% Base signal, start with the discharge pulse. Discharge assumed as
% positive
if cDmax >= cCmax
    baseSig = [C1*ones(round(T1*fs),1); zeros(round(T2*fs),1); -C2*ones(round(T3*fs),1); zeros(round(T4*fs),1)];
else
    baseSig = [C2*ones(round(T3*fs),1); zeros(round(T2*fs),1); -C1*ones(round(T1*fs),1); zeros(round(T4*fs),1)];
end

% Number of samples per period
N  = length(baseSig);

% FFT of base-signal
BS = fft(baseSig);
harmonics = [0:N-1];        % General harmonic set

% Find non excited and excited base-signal harmonics
idx_supp = db(BS) < -100;
idx_exc = db(BS) > -100;

suppBaseSigHarm = harmonics(idx_supp); % Suppressed base-signal harmonics
excBaseSigHarm = harmonics(idx_exc);   % Excited base-signal harmonics

Hlb = 1;                    % Minimum harmonic
Hub = floor(N*fMax/fs);     % Maximum harmonic corresponding to fMax

harmSet = [Hlb:Hub]';        % General harmonics set within bandwidth
harmSupp = intersect(harmSet,suppBaseSigHarm); % Suppressed harmoncis within bandwidth
harmExc = intersect(harmSet,excBaseSigHarm);   % Excited harmoncis within bandwidth

% Generate random phase multisine
K = beta*cMin;                          % Multisine scale factor
msSignal = MultisineSignal(N,harmExc,K);

% Pulse-multisine and time vector, a structure variable
u.pmSignal = cRate*(baseSig+msSignal);
u.timeVec = [0:N-1]'/fs;
% Propeties of the pulse-multisine
u.harmSupp = harmSupp;
u.harmExc = harmExc;
u.cDmax = cDmax;
u.cCmax = cCmax;
u.T1 = T1;
u.T3 = T3;
u.T2 = T2;
u.T4 = T4;
u.N = N;
u.alpha = alpha;
u.fMax = fMax;
u.fs = fs;
u.cRate = cRate;
u.baseSignal = baseSig*cRate;
u.msSignal = msSignal*cRate;

end

function msSignal = MultisineSignal(N, harmExc, K)
% Generate a random phase multisine

% Number of excited harmonics
lU = length(harmExc);

% Random phases
randPhase = 2*pi*rand(lU,1);

%Frequency components
Au = exp(1i*randPhase);

linesExc = harmExc + 1;

% Initialise multisine FFT
msSignalFFT = zeros(N,1); 

% Populate excited components
msSignalFFT(linesExc) = Au;

% Time signals
msSignal = 2*real(ifft(msSignalFFT));

% Scale multisine signal
msSignal = msSignal/max(abs(msSignal))*K;

end
