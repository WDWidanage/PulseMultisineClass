classdef PMObj < handle
    % Class to design pulse-multisine signals and estimate a Nonlinear ECM based on pulse-multisine data
    %
    % References:
    %  Widanage, W. D., Barai, A., Chouchelamane, G.H., Uddin, K., McGordon, A.,
    %  Marco, J. and Jennings, P., "Design and use of multisine signals for Li-ion battery equivalent circuit modelling. Part 1: Signal design",
    %  Journal of Power Sourcers, 324, pp. 70-78.
    %
    %  Widanage, W. D., Barai, A., Chouchelamane, G.H., Uddin, K., McGordon, A.,
    %  Marco, J. and Jennings, P., "Design and use of multisine signals for Li-ion battery equivalent circuit modelling. Part 2: Model estimation",
    %  Journal of Power Sourcers, 324, pp. 61-69.
    %
    % Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 01/07/2016 (For whom the bells toll!)
    % All Rights Reserved
    % Software may be used freely for non-comercial purposes only
    
    properties
        % General properties of reference cell, reference SoC and
        % temperature at which a pulse-multisine is genearated or the
        % measured at.
        refCell = [];              % Cell manufacture and details. Input should be of type string
        refSoC = [];               % SoC at which pulse-multisine is expected to be applied or collected
        refTemp = [];              % Temperature at which pulse-multisine is expected to be applied or collected
        
        warningFlag = [];          % If any warnings arise a flag will be generated. This can take a value between 1-3 indicating
        %  1 - Flag high transients in voltage. Voltage may not have reached a sufficient steady-state level. LPM transient option should be set to 1. 'obj'.estSetting.LPMTransient = 1.
        %  2 - Flag if model poles are unstable or complex. Model order or frequency range should be reduced via 'obj'.estSetting.ECMorder or 'obj'.estSetting.fCutOff.
        %  3 - Flag rank deficiency in regressor matrix. Model order or frequency range should be reduced via 'obj'.estSetting.ECMorder or 'obj'.estSetting.fCutOff.
        
        % Properties of pulse-multisine signal parameters and null feilds for the generated signal
        refSig = struct(...
            'cRate', [],...        % C-rate of the battery
            'cDmax', [],...        % Maximum applicable 10 s discharge current, for a given SoC and temperature, as specified by the manufacture
            'cCmax', [],...        % Maximum applicable 10 s dharge current, for a given SoC and temperature, as specified by the manufacture
            'T1',    [],...        % Time interval of the largest base-signal pulse
            'T2',    [],...        % Time interval of the first base-signal rest period
            'T4',    [],...        % Time interval of second base-signal rest period
            'alpha', [],...        % Fraction of smaller base-signal pulse compared to maximum allowed c-rate
            'fMax',  [],...        % Highest excited frequency in the multisine signal. This should be set to cover the frequency of interest. 1Hz is sufficient for a drive-cycle
            'fs',    [],...        % Sampling frequency at which the cell cylcer is run. Normmaly it is at 10Hz. Note that fs should always be fs >= 2*fMax.
            'pmSignal',   [],...   % Generated pulse-multisine will be saved in this
            'timeVec',    [],...   % Signal time column
            'harmSupp',   [],...   % List of any suppressed harmonics
            'harmExc',    [],...   % List of excited harmonics
            'T3',         [],...   % Duration of largest base-signal pulse
            'baseSignal', [],...   % Generated base-signal
            'msSignal',   [],...   % Multisine signal
            'N',          []);     % Signal period length in samples
        
        % Properties of measured voltage and current
        measSig = struct(...
            'measCurr', [],...    % place holder for measured current (A)
            'measVol',  [],...    % place holder for measured voltage (V)
            'measTime', [],...    % place holder for measured time (s)
            'P',[]);              % Number of measured periods
        
        % Properties of parameter estimation settings
        estSetting = struct(...
            'LPMTransient',  0,...          % Specify if tranasient term in impedance should be estimated
            'LPMi',         'n',...         % Use LPMi if input is a non-zero mean signal
            'fCutOff',       1,...
            'useElis',      'n',...         % Use FDIDENT elis for ECM parameterisation
            'JacobianECM',  'on',...
            'JacobianSig',  'on',...
            'ECMOrder',      2,...          % Number of RC pairs
            'startPeriod',   2);            % Start from this period for FRF estimation
        
        % Place holder for estimated battery impedance
        estimatedImpedance = struct(...
            'freq',       [],...           % Frequency list
            'impedance',  [],...           % Estimated impedance
            'impVar',     []);             % Variance of impedance
        
        % Place holder for NLECM parameters
        nlECMPara = struct(...
            'Ro',          [],...   % ECM internal resistance
            'Rp',          [],...   % ECM polarisation resistances
            'Tau',         [],...   % ECM time constants
            'Cp',          [],...   % ECM polarisation capacitances
            'c1',          [],...   % Nonlinear sigmoid coefficient c1
            'c2',          [],...   % Nonlinear sigmoid coefficient c2
            'ecmStd',      [],...   % ECM parameter standard deviations
            'nlStd',       [],...   % Non-linear sigmoid parameter standard deviations
            'cfECM',       [],...   % ECM cost-function at termination
            'cfSigmoid',   [],...   % Non-linear sigmoid cost-function at termination
            'ecmAIC',      []);     % ECM AIC value
    end
    
    
    methods
        % Constructor method
        function obj = PMObj(varargin)
            
            % Constructor method for pulse-multisine design
            if nargin == 1
                obj.refSig = varargin{1};
                if isfield(varargin{1},'refSoC')
                    obj.refSoC = varargin{1}.refSoC;
                end
                if isfield(varargin{1},'refTemp')
                    obj.refTemp = varargin{1}.refTemp;
                end
                if isfield(varargin{1},'refCell')
                    obj.refCell = varargin{1}.refCell;
                end
                
                obj = obj.createPMSignal();
                
            end
            
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Pulse-multisine generating method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = createPMSignal(obj)
            u = PulseMultisine(obj.refSig.cDmax, obj.refSig.cCmax, obj.refSig.cRate, obj.refSig.fs, obj.refSig.T1, obj.refSig.T2, obj.refSig.T4, obj.refSig.alpha, obj.refSig.fMax);
            obj.refSig.pmSignal = u.pmSignal;
            obj.refSig.timeVec = u.timeVec;
            obj.refSig.harmSupp = u.harmSupp;
            obj.refSig.harmExc = u.harmExc;
            obj.refSig.T3 = u.T3;
            obj.refSig.baseSignal = u.baseSignal;
            obj.refSig.msSignal = u.msSignal;
            obj.refSig.N = u.N;
        end
        
        % funtion set a pre-designed pulse-multisine data strucutre
        function obj = setRefSig(obj, u)
            refFieldNames = fieldnames(obj.refSig);
            uFieldNames = fieldnames(u);
            
            for ii = 1:length(refFieldNames)
                [~,idx] = ismember(refFieldNames{ii},uFieldNames);
                if ~isempty(idx)
                    obj.refSig.(refFieldNames{ii}) = u.(uFieldNames{idx})(:);
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameter estiamtion methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = estNonParaImpedance(obj)
            % Non-parametric impedance estimation
            startIdx = 1;
            if isempty(obj.refSig.N)
                N = length(obj.measSig.measCurr);
            else
                N = obj.refSig.N;
            end
            
            endIdx = N * obj.measSig.P;
            dataSeg = [startIdx:endIdx];
            currVec = obj.measSig.measCurr(dataSeg);
            volVec = obj.measSig.measVol(dataSeg);
            fs = obj.refSig.fs;
            linesExc = obj.refSig.harmExc + 1;
            freqExc = obj.refSig.harmExc/N*fs;
            
            P = obj.measSig.P;
            
            if P > 1
                startPeriod = obj.estSetting.startPeriod;
                optionsLPM.transient = obj.estSetting.LPMTransient;
                
                currAllPeriods = reshape(currVec,N,P);
                volAllPeriods = reshape(volVec,N,P);
                
                % Check level of transient
                lastPeriodRms = (rms(volAllPeriods(:,P)))^2;
                firstPeriodRms = (rms(volAllPeriods(:,startPeriod)))^2;
                perVarDiff = (firstPeriodRms - lastPeriodRms)/lastPeriodRms*100;
                
                
                if perVarDiff > 20 && ~obj.estSetting.LPMTransient
                    msg = sprintf('High transient error detected. Consider either setting: \n Starting period larger than %g. ''obj''.estSetting.startperiod \n Setting LPM tranient to 1. ''obj''.estSetting.LPMTransient = 1\n',startPeriod);
                    warning(msg);                           % Flag high transients in voltage.
                    obj.warningFlag(1,1) = 1;
                else
                    obj.warningFlag(1,1) = 0;
                end
                
                % Calculate mean over periods
                currMean = mean(currAllPeriods(:,startPeriod:P),2);
                volMean = mean(volAllPeriods(:,startPeriod:P),2);
                
                volMean0  = volMean - mean(volMean);
                
                % Estimate Impedance
                
                
                % FFT
                U = fft(currMean);
                Y = fft(volMean0);
                
                % call LPM
                [Glpm, ~, ~ ,CGlpm] = LPM(U(linesExc),Y(linesExc),linesExc,optionsLPM);
                
                obj.estimatedImpedance.freq = freqExc;
                obj.estimatedImpedance.impedance = Glpm;
                obj.estimatedImpedance.impVar = CGlpm;
            else                                             % Non-periodic measured data
                if ismember(lower(obj.estSetting.LPMi),'n')  % Use LPM with transients
                    
                    optionsLPM.transient = obj.estSetting.LPMTransient;
                    
                    % FFT
                    U = fft(currVec);
                    Y = fft(volVec);
                    
                    % call LPM
                    [Glpm, ~, ~ ,CGlpm] = LPM(U(linesExc),Y(linesExc),linesExc,optionsLPM);
                    
                    obj.estimatedImpedance.freq = freqExc;
                    obj.estimatedImpedance.impedance = Glpm;
                    obj.estimatedImpedance.impVar = CGlpm;
                    
                else                                         % Use LPMi with transients
                    est = LPMi(currVec,volVec,linesExc,fs);
                    obj.estimatedImpedance.freq = freqExc;
                    obj.estimatedImpedance.impedance = est.G;
                    obj.estimatedImpedance.impVar = est.Cg;
                end
            end
            
            
            
        end
        
        
        
        function [obj,freqExc_bw, G_bw, GTF] = estImpedance(obj)
            % Impedance parameterisation
            
            % First call impedance estimation method
            obj = obj.estNonParaImpedance();
            
            % Start of impedance parameterisation
            freqExc = obj.estimatedImpedance.freq;
            fCutOff = obj.estSetting.fCutOff;
            G = obj.estimatedImpedance.impedance;
            V = obj.estimatedImpedance.impVar;
            nb = obj.estSetting.ECMOrder;
            fs = obj.refSig.fs;
            
            
            if isnan(fCutOff)
                idx_bw = freqExc <= 1;
            elseif isempty(fCutOff)
                idx_bw = freqExc <= 1;
            else
                idx_bw = freqExc <= fCutOff;
            end
            
            
            G_bw = G(idx_bw);
            varG_bw = V(idx_bw);
            freqExc_bw = freqExc(idx_bw);
            
            na = nb;
            wNorm = 2*pi*freqExc_bw'/fs;
            
            
            if strcmpi(obj.estSetting.useElis,'y')
                fData = fiddata(G_bw,ones(size(G_bw)),freqExc_bw,varG_bw);
                itMax = 100;
                G_fit = elis(fData,'s',nb,na,struct('algorithm','LM','fs',fs,'itmax',itMax,'stabilization','','forceminimumphase','')); % stabilisation and forceminimumphase argument is 'r'
                %Extract coefficients
                Bn =  G_fit.num./G_fit.denom(end);
                An =  G_fit.denom./G_fit.denom(end);
                
                [residueTF,polesTF,RoTF] = residue(Bn,An); % Partial fraction expansion
                
                
                % Store ECM parameters,
                obj.nlECMPara.Ro = -RoTF;
                obj.nlECMPara.Rp = residueTF./polesTF;
                obj.nlECMPara.Tau = -1./polesTF;
                obj.nlECMPara.Cp = obj.nlECMPara.Tau ./obj.nlECMPara.Rp;
                
                % Model ECM frf
                opt.fs = fs;
                Gsplit = ECMFRF([obj.nlECMPara.Ro,obj.nlECMPara.Rp', obj.nlECMPara.Tau'],wNorm,opt);
                lW = length(wNorm);
                GTF = Gsplit(1:lW)+1i*Gsplit(lW+1:end);
                
                
            else
                
                % Initial TF estimation using Levi method
                leviOptions.fs = fs;
                Model = LeviAlgorithm(G_bw,wNorm,nb,na,leviOptions); % Initial TF parameter guess
                
                
                % Map the Levi TF parameters to the ECM parameters
                [residueTF,polesTF,RoTF] = residue(Model.B,[Model.A;1]);  % Partial fraction expansion
                
                if isempty(RoTF)
                    RoTF = 0;
                end
                RpTF = -residueTF./polesTF;
                TauTF = -1./polesTF;
                
                % ECM parameters from Levi method as intial values for
                RpTF = -(abs(RpTF)); % Due to numerical isues in the residue function, force RpTF values to be all negative. Like Ro, Rp should be negative due to the negative transfer function gain (discharge current is assumed positve while voltage is decreasing)
                thetaECM0 = real([RoTF,RpTF',TauTF']);                      % Assume only the resistive parts for intial guess for non-linear optimisation
                
                % Optimisation of ECM parameters via LM method
                opt.fs = fs;
                fhECMFRF = @(theta,w)ECMFRF(theta,w,opt);
                G_bw_real = real(G_bw);
                G_bw_imag = imag(G_bw);
                Gd = [G_bw_real;G_bw_imag];
                
                optionsECM.s = sqrt([varG_bw; varG_bw]/2);
                optionsECM.Jacobian = obj.estSetting.JacobianECM;
                [tfECMParaOpt,infoECMTFPara] = LMAlgorithm(fhECMFRF, thetaECM0, wNorm, Gd, optionsECM);
                
                % Model ECM frf
                Gsplit = ECMFRF(tfECMParaOpt,wNorm,opt);
                lW = length(wNorm);
                GTF = Gsplit(1:lW)+1i*Gsplit(lW+1:end);
                
                % Corrected AIC
                dof = length(tfECMParaOpt)+1;
                dtPts = length(Gd);
                
                
                % Store ECM parameters, AIC and cost-fucntion value
                obj.nlECMPara.Ro = -tfECMParaOpt(1);
                [obj.nlECMPara.Rp,idx] = sort(-tfECMParaOpt(2:nb+1));
                Tau =  tfECMParaOpt(nb+2:end);
                obj.nlECMPara.Tau = Tau(idx);
                obj.nlECMPara.Cp = obj.nlECMPara.Tau./obj.nlECMPara.Rp;
                RoStd = infoECMTFPara.stdTheta(1);
                RpStd = infoECMTFPara.stdTheta(2:nb+1);
                TauStd = infoECMTFPara.stdTheta(nb+2:end);
                obj.nlECMPara.ecmStd = [RoStd; RpStd(idx);TauStd(idx)];     
                obj.nlECMPara.ecmAIC = dtPts*log(infoECMTFPara.cF_iter(end)/dtPts) + 2*dof + 2*dof*(dof+1)/(dtPts-dof-1);
                obj.nlECMPara.cfECM = infoECMTFPara.cF_iter(end);
                
                % Calculate ECM poles
                polesECM = -1./obj.nlECMPara.Tau;
                
                
                % Generate warining cases
                msgFix1 = sprintf('\n Consider using a lower ECM model order than %g. ''obj''.estSetting.ECMOrder', obj.estSetting.ECMOrder);
                msgFix2 = sprintf('\n Consider reducing the fitting frequency range to less than %g. ''obj''.estSetting.fCutoff \n',obj.refSig.fMax);
                
                if any(polesECM  > 0) && obj.estSetting.ECMOrder > 1                 % Flag warning if poles are unstable
                    warning(['Unstables ECM poles',msgFix1,msgFix2]);
                    obj.warningFlag(2,1) = 2;
                elseif  any(polesECM  > 0)
                    warning(['Unstables ECM poles',msgFix2]);
                    obj.warningFlag(2,1) = 2;
                else
                    obj.warningFlag(2,1) = 0;
                end
                
                if ~isreal(polesECM) && obj.estSetting.ECMOrder > 1                  % Flag warning if poles are complex
                    warning(['Complex  ECM poles',msgFix1,msgFix2]);
                    obj.warningFlag(2,1) = 2;
                elseif  ~isreal(polesECM)
                    warning(['Complex ECM poles',msgFix2]);
                    obj.warningFlag(2,1) = 2;
                else
                    obj.warningFlag(2,1) = 0;
                end
                
                
                if any(~infoECMTFPara.LMRankFull)  && obj.estSetting.ECMOrder > 1    % Flag LM rank deiciency
                    warning(['LM Regressor rank deficient',msgFix1,msgFix2]);
                    obj.warningFlag(3,1) = 3;
                elseif any(~infoECMTFPara.LMRankFull)
                    warning(['LM Regressor rank deficient',msgFix2]);
                    obj.warningFlag(3,1) = 3;
                else
                    obj.warningFlag(3,1) = 0;
                end
            end
        end
        
        
        function [obj,yECM0,volMean0,yNL0,volMean,yNL] = estNLCharac(obj)
            % Non-linear characterisation
            
            startIdx = 1;
            endIdx = obj.refSig.N * obj.measSig.P;
            dataSeg = [startIdx:endIdx];
            timeVec = obj.measSig.measTime(dataSeg);
            currVec = obj.measSig.measCurr(dataSeg);
            volVec = obj.measSig.measVol(dataSeg);
            startPeriod = obj.estSetting.startPeriod;
            freqExc = obj.estimatedImpedance.freq;
            linesExc = obj.refSig.harmExc + 1;
            nb = obj.estSetting.ECMOrder;
            N = obj.refSig.N;
            P = obj.measSig.P;
            fs = obj.refSig.fs;
            fCutOff = obj.estSetting.fCutOff;
            
            
            
            
            currAllPeriods = reshape(currVec,N,P);
            volAllPeriods = reshape(volVec,N,P);
            
            % Calculate mean over periods
            currMean = mean(currAllPeriods(:,startPeriod:P),2);
            volMean = mean(volAllPeriods(:,startPeriod:P),2);
            OCV = mean(volMean);
            volMean0  = volMean - OCV;
            
            % FFT
            U = fft(currMean);
            
            if isnan(fCutOff)
                idx_bw = freqExc <= 1;
            elseif isempty(fCutOff)
                idx_bw = freqExc <= 1;
            else
                idx_bw = freqExc <= fCutOff;
            end
            
            
            freqExc_bw = freqExc(idx_bw);
            wNorm = 2*pi*freqExc_bw/fs;
            linesExc_bw = linesExc(idx_bw);
            lW = length(wNorm);
            
            % Estimate steady state linear over-voltage signal
            paraVec = [-obj.nlECMPara.Ro;-obj.nlECMPara.Rp;obj.nlECMPara.Tau];
            opt.fs = fs;
            Gsplit = ECMFRF(paraVec,wNorm,opt);
            GECM = Gsplit(1:lW)+1i*Gsplit(lW+1:end);
            
            YECM0  = (GECM.*U(linesExc_bw));
            
            yECM0Tmp = zeros(N,1);
            yECM0Tmp(linesExc_bw) = YECM0;
            yECM0 = 2*real(ifft(yECM0Tmp));
            
            
            % Esimate Sigmoid coefficients
            optionsSig.Jacobian = obj.estSetting.JacobianSig;
            [sigParaOpt,infoSigPara] = LMAlgorithm(@SigmoidFcn,[0.5,0.1],yECM0,volMean0,optionsSig);
            yNL0 = SigmoidFcn(sigParaOpt,yECM0);
            
            
            % Simulate overall model
            thetaECM = [obj.nlECMPara.Ro; obj.nlECMPara.Rp; obj.nlECMPara.Tau];
            volLin = ECMStable(thetaECM,currVec,0,timeVec,nb);
            volNL = SigmoidFcn(sigParaOpt,volLin) + OCV;
            
            % Eliminate first period and average
            simVolAllPeriods = reshape(volNL,N,P);
            yNL = mean(simVolAllPeriods(:,startPeriod:P),2);
            
            
            obj.nlECMPara.c1 = sigParaOpt(1);
            obj.nlECMPara.c2 = sigParaOpt(2);
            obj.nlECMPara.nlStd = infoSigPara.stdTheta;
            obj.nlECMPara.cfSigmoid = infoSigPara.cF_iter(end);
            
        end
        
        
        function obj = estNLECM(obj)
            % Parameterise full NL-ECM
            obj = obj.estImpedance;
            obj = obj.estNLCharac;
        end
        
        %%%%%%%%%%%%%%%%%
        % Ploting methods
        %%%%%%%%%%%%%%%%%
        function obj = plotRefSig(obj)
            cRate = obj.refSig.cRate;
            fs = obj.refSig.fs;
            maxIrate = max(obj.refSig.pmSignal)/cRate;
            minIrate = min(obj.refSig.pmSignal)/cRate;
            currRms = rms(obj.refSig.pmSignal);
            BS = fft(obj.refSig.baseSignal);
            
            freqRec = [0:obj.refSig.N-1]/(obj.refSig.N)*fs;
            
            U = fft(obj.refSig.pmSignal);
            excitedLines = (obj.refSig.harmExc)+1;        % Excited lines realtive to record
            freqExc = (excitedLines-1)*fs/(obj.refSig.N);
            freqSupp = obj.refSig.harmSupp*fs/obj.refSig.N;
            
            
            figure()
            plot(obj.refSig.timeVec,obj.refSig.pmSignal,'- .');
            xlabel('Time (s)'); ylabel ('Current (A)'); title(['Cmin: ',num2str(minIrate),' Cmax: ',num2str(maxIrate), ' Rms: ',num2str(currRms)])
            
            figure()
            hist(obj.refSig.pmSignal,18)
            xlabel('Signal amplitude'); ylabel('Frequency of occurrence')
            
            % Perform Coulomb counting to check SoC variation
            soc= CoulombCounting(obj.refSig.pmSignal,obj.refSig.timeVec,cRate,cRate,obj.refSoC/100);
            socChange = [max(soc)-min(soc)]*100;
            
            figure()
            plot(obj.refSig.timeVec,soc*100)
            xlabel('Time (s)'); ylabel('SoC (%)'); title(['SoC Change: ', num2str(socChange),'%'])
            
            figure();
            plot(freqRec,abs(BS),'-x',freqExc,abs(U(excitedLines)),'-or',freqSupp,abs(U(obj.refSig.harmSupp+1)),'go');
            xlabel('Frequency (Hz)')
            ylabel('FFT magnitude (abs)')
            xlim([0, 1])
            
        end
        
        function obj = plotMeasSig(obj)
            startIdx = 1;
            if isempty(obj.refSig.N)
                N = length(obj.measSig.measCurr);
            else
                N = obj.refSig.N;
            end
            endIdx = N * obj.measSig.P;
            dataSeg = [startIdx:endIdx];
            
            timeVec = obj.measSig.measTime(dataSeg) - obj.measSig.measTime(startIdx); % Set time to start from zero
            currVec = obj.measSig.measCurr(dataSeg);
            volVec = obj.measSig.measVol(dataSeg);
            
            % plot measured and reference current signal
            refCurr = repmat(obj.refSig.pmSignal,obj.measSig.P,1);
            if isempty(refCurr)
                refCurr = zeros(size(currVec));
                refMeasCurrErr = zeros(size(currVec));
            else
                refMeasCurrErr = currVec - refCurr;
            end
            
            if obj.measSig.P > 1
                volAllPeriods = reshape(volVec,obj.refSig.N,obj.measSig.P);
                volTransAllPeriods = volAllPeriods - repmat(volAllPeriods(:,obj.measSig.P),1,obj.measSig.P);
                volTrans = volTransAllPeriods(:);
                
                % Plot measured signals
                figure()
                subplot(4,1,1)
                plot(timeVec,currVec,timeVec,refCurr,'.')
                xlabel('Time (s)'); ylabel('Current (A)'); legend('Measured','Reference')
                subplot(4,1,2)
                plot(timeVec,volVec)
                xlabel('Time (s)'); ylabel('Voltage (V)');
                subplot(4,1,3)
                plot(timeVec,refMeasCurrErr)
                xlabel('Time (s)'); ylabel('Current error (A)');
                subplot(4,1,4)
                plot(timeVec,volTrans)
                xlabel('Time (s)'); ylabel('Voltage transients (V)');
            else
                figure()
                subplot(3,1,1)
                plot(timeVec,currVec,timeVec,refCurr,'.')
                xlabel('Time (s)'); ylabel('Current (A)'); legend('Measured','Reference')
                subplot(3,1,2)
                plot(timeVec,volVec)
                xlabel('Time (s)'); ylabel('Voltage (V)');
                subplot(3,1,3)
                plot(timeVec,refMeasCurrErr)
                xlabel('Time (s)'); ylabel('Current error (A)');
            end
            
        end
        
        % Plot impedance
        function plotImpedance(obj)
            G = obj.estimatedImpedance.impedance;
            f = obj.estimatedImpedance.freq;
            V = obj.estimatedImpedance.impVar;
            
            figure()
            subplot(2,1,1)
            plot(f,db(G),f,db(V)/2)
            legend('Impedance','Variance')
            xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)')
            subplot(2,1,2)
            plot(f,180/pi*unwrap(angle(G)))
            xlabel('Frequency (Hz)'); ylabel('Phase (deg)')
        end
        
        
        % Plot estimated impedance and fitted ECM
        function plotECMFit(obj)
            [obj, freqExc_bw, G_bw, GTF] = obj.estImpedance;
            
            figure()
            subplot(2,1,1)
            plot(freqExc_bw,db(G_bw),'o',freqExc_bw,db(GTF),'. -')
            legend('Estimated impedance','Fit')
            xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)')
            titleTF = sprintf(['Cost function %E. AIC: %g'],obj.nlECMPara.cfECM,obj.nlECMPara.ecmAIC);
            title(titleTF)
            subplot(2,1,2)
            plot(freqExc_bw,180/pi*unwrap(angle(G_bw)),'o',freqExc_bw, 180/pi*unwrap(angle(GTF)),'. -')
            xlabel('Frequency (Hz)'); ylabel('Phase (deg)')
        end
        
        function plotNLfit(obj)
            [obj, yECM0,volMean0,yNL0] = obj.estNLCharac;
            
            figure()
            [~, idxSort] = sort(yECM0);
            plot(yECM0(idxSort),volMean0(idxSort),'r o',yECM0(idxSort),yNL0(idxSort));
            xlabel('Model over-voltage (V)'); ylabel('Measured over-voltage (V)')
            titleSig = sprintf(['Linear coeff: %E   NL coeff: %E  Cost Fcn: %E'],obj.nlECMPara.c1,obj.nlECMPara.c2,obj.nlECMPara.cfSigmoid);
            title(titleSig)
        end
        
        function obj = plotSimVoltage(obj)
            % Plot simulated NLECM voltage with measured voltage
            [obj,~,~,~,volMean,yNL] = estNLCharac(obj);
            timeOnePeriod = [0:obj.refSig.N-1]*obj.refSig.fs;
            
            % Calculate R^2 value. Coeffcicnet of determination
            ssRes = var(volMean  - yNL);
            ssTot = var(volMean);
            rSq = (1 - ssRes/ssTot)*100;
            
            figure()
            plot(timeOnePeriod,volMean,'o',timeOnePeriod,yNL,'. -')
            xlabel('Time (s)'); ylabel('Voltage (A)'); legend('Measured','Simulated')
            title(sprintf('R squared: %2.1f',rSq))
        end
        
        function plotAll(obj)
            % Plot all the figures from each stage of analysis
            obj.plotMeasSig;
            obj.plotECMFit;
            obj.plotNLfit;
            obj.plotSimVoltage;
        end
    end
    
end

