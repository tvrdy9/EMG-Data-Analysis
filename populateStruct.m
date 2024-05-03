clearvars -except VIBdata;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code reads in all force/emg files to populate structure VIBdata %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% load data structure and CV indeces from Excel file %%%
fsamp = 2000;
f=1;    % file index in excel sheet
CVStartIndex = uint32(readmatrix('E:/Vibration/Participant_data.xlsx','Sheet','General','Range','J2:J169'));
TrialNum = double(readmatrix('E:/Vibration/Participant_data.xlsx','Sheet','General','Range','D2:D169'));
epochSize = 1;
stepSize = epochSize/4;
numEpochs = (10-epochSize)/stepSize + 1;

% %1-2 Hz alignment - run first time for manual alignment
% fc = ones(168,2) + [zeros(168,1),ones(168,1)];
% forceLPF = 2*ones(168,1);
% run = 1;

%1-2 Hz alignment - run 2nd time to save EMdelays and optimal filter frequencies to excel
EMdelaySamples = double(readmatrix('E:/Vibration/Participant_data.xlsx','Sheet','General','Range','K2:K169'));
fc = ones(168,2) + [zeros(168,1),ones(168,1)];
forceLPF = 2*ones(168,1);
run = 2;

% %0.75 - 5 Hz alignment - run 3rd time to populate struct with all vars
% EMdelaySamples = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','G2:G169'));
% fc = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','E2:F169'));
% forceLPF = 5*ones(168,1);
% run = 3;

% %fcHPF - fcLPF Hz alignment - run 3rd time to populate struct with all vars
% EMdelaySamples = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','AA2:AA169'));
% fc = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','Y2:Z169'));
% forceLPF = fc(:,2);
% run = 3;

% %1.5 - 5 Hz alignment - run 3rd time to populate struct with all vars
% EMdelaySamples = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','AU2:AU169'));
% fc = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','AS2:AT169'));
% forceLPF = 5*ones(168,1);
% run = 3;

%%% loop through all subjects %%%
for s = 1:length(VIBdata.subject)
    %%% loop through each condition %%%
    for c = 1:3
        switch c
            case 1
                conditionString = "Vision";
            case 2
                conditionString = "Sound";
                if TrialNum(f) == 0 %check for vision only subjects
                    f = f+4;
                    break;
                end
            case 3
                conditionString = "Vibration";
                if TrialNum(f) == 0 %check for vision only subjects
                    f = f+2;
                    break;
                end
        end
        %%% populate individual 5 and 20 percent %%%
        for t=[5,20]
            
            % load EMG file
            filename = sprintf('%s-%s-%d-Final.mat',VIBdata.subject(s).data.id,conditionString,t);
            pathname = sprintf('E:/Vibration/Final/%s',VIBdata.subject(s).data.id);
            filename = fullfile(pathname,filename);
            load(filename);     % load variables into workspace
            
            % load Force file
            filename = sprintf('%s-%s-%d-Force.mat',VIBdata.subject(s).data.id,conditionString,t);
            pathname = sprintf('E:/Vibration/Force/%s',VIBdata.subject(s).data.id);
            filename = fullfile(pathname,filename);
            tempForce = load(filename); % load variables into struct
            
            % calculate force variables
            force = tempForce.force_save(1:SIGlength*fsamp); % save force values after trigger
            steadyForceInd = [CVStartIndex(f) CVStartIndex(f)+10*fsamp]; % steady force indeces
            
            % calculate emg variables
            steadyEMGInd = [CVStartIndex(f) CVStartIndex(f)+10*fsamp]; %steady emg indeces
            numSamples = floor(fsamp*SIGlength);
            CST = zeros(1,numSamples);
            binary = zeros(size(MUPulses,2),numSamples);

            % populate CST
            for mu = 1:numel(MUPulses)
                for spikeIndex = 1:numel(MUPulses{mu})
                    tempSpikeTime = MUPulses{mu}(spikeIndex);
                    binary(mu,tempSpikeTime) = 1;
                    CST(tempSpikeTime) = CST(tempSpikeTime) + 1;
                end
            end

            % calculate active index and normalize CST
            win = 0.25*fsamp; % 250 ms window
            activeind = zeros(size(binary,1),size(binary,2));
            for r = 1:size(binary,1)
                start = 1;
                endd = size(binary,2);
                num = floor((endd-start)/win)-1;
                for i = 0:num-1 % slide 250 ms window to calculate active index
                    ind = start+(win*i);
                    temp = binary(r,ind:ind+win-1);
                    if sum(temp) > 0
                        activeind(r,ind:ind+win-1) = 1;
                    elseif sum(temp) == 0
                        activeind(r,ind:ind+win-1) = 0;
                    end
                    if i==num-1
                        ind = start+(win*(i+1));
                        temp = binary(r,ind:end);
                        if sum(temp) > 0
                            activeind(r,ind:end) = 1;
                        elseif sum(temp) == 0
                            activeind(r,ind:end) = 0;
                        end
                    end
                end
            end

            numMU = sum(activeind);
            activeind = numMU;
            normCST = CST./activeind;
            normCST(isinf(normCST)|isnan(normCST)) = 0;
            
            % calculate ISIs
            ISIs = cell(1,size(MUPulses,2));
            for mu = 1:numel(MUPulses)
                tempISI = zeros(1,numel(MUPulses{mu})-1);
                for spikeIndex = 2:numel(MUPulses{mu})
                    tempISI(1,spikeIndex-1) = MUPulses{mu}(spikeIndex)-MUPulses{mu}(spikeIndex-1);
                end
                ISIs{mu} = tempISI;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% calculate variables over steady 10 seconds              %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            numMUs = numel(MUPulses);
            numSteadyEMGSamps = 10*fsamp;
            CVforce = zeros(1,numEpochs);
            SDforce = zeros(1,numEpochs);
            CVISI = zeros(numMUs,numEpochs); 
            SDCST = zeros(1,numEpochs);
            usableMUs = zeros(1,numEpochs);
            steadyPNRs = zeros(1,numMUs);
            
            %smooth CST and filter force
            hann_window = hann(floor(0.4*fsamp));
            smoothCST = conv(normCST,hann_window,'same');
            fcForce = 20;
            [b20, a20] = butter(6,fcForce/fsamp/2);
            force = filtfilt(b20,a20,force);
            
            %%% correct for EMdelay between EMG and Force signals %%%
            fcHPF = fc(f,1);
            fcLPF = fc(f,2);
            [bHPF,aHPF] = butter(2,fcHPF/(fsamp/2),'high');
            [bLPF,aLPF] = butter(2,fcLPF/(fsamp/2));
            [bLPFforce,aLPFforce] = butter(2,forceLPF(f)/(fsamp/2));
            filtCST = filtfilt(bLPF,aLPF,smoothCST);
            filtCST = filtfilt(bHPF,aHPF,filtCST);
            filtForce = filtfilt(bHPF,aHPF,force);
            filtForce = filtfilt(bLPFforce,aLPFforce,filtForce);
            iEMG = double(SIG{8,4})';
            iEMG = bandpass(iEMG,[20 500],fsamp);
            iEMG = abs(iEMG);
            [bEMG,aEMG] = butter(2,2/(fsamp/2));
            iEMG = filtfilt(bEMG,aEMG,iEMG);
            
            if (run == 1)
                %manually find EM delay and cross-correlate for 100 ms fine tuning
                EMdelay = -1*round(delaySlider(fsamp, iEMG, force, filtCST, filtForce, steadyEMGInd, steadyForceInd, titleString));
                steadyEMGInd = steadyEMGInd - EMdelay;
                maxDelay = 0.05*fsamp; %max of 100 ms correction
                fluctForce = filtForce(steadyForceInd(1):steadyForceInd(2)-1);
                fluctCST = filtCST(steadyEMGInd(1):steadyEMGInd(2)-1);
                fineDelay = finddelay(fluctForce,fluctCST,maxDelay);
                EMdelay = EMdelay - fineDelay;
                steadyEMGInd = steadyEMGInd + fineDelay;
%                 cellString = sprintf('K%d',f+1);
%                 filename = "Participant_data.xlsx";
%                 pathname = "E:/Vibration";
%                 filename = fullfile(pathname,filename);
%                 writematrix(EMdelay,filename,'Sheet','General','Range',cellString);
                
            elseif (run == 2)
                %%% populate excel for [0.75 - 5 Hz] force bandwidth %%%
                %read in EMdelays from excel
                EMdelay = EMdelaySamples(f);
                tempSteadyEMGInd = steadyEMGInd - EMdelay;
                
                %%% find optimal cutoff freqs %%%
                maxDelay = 0.2*fsamp; %max of 200 ms correction
                
                %find max cross-correlation indeces for 0.75 - 5 Hz bandwidth
%                 fcHPF = 0.75;
                fcHPF = 1.5;
                fcLPF = 5;
                [bHPF,aHPF] = butter(2,fcHPF/(fsamp/2),'high');
                HPFfiltForce = filtfilt(bHPF,aHPF,force);
                HPFfiltCST = filtfilt(bHPF,aHPF,smoothCST);
                [bLPF,aLPF] = butter(2,fcLPF/(fsamp/2));
                filtForce = filtfilt(bLPF,aLPF,HPFfiltForce);
                fluctForce = filtForce(steadyForceInd(1):steadyForceInd(2)-1);
                normFluctForce = fluctForce/max(abs(fluctForce));
                filtCST = filtfilt(bLPF,aLPF,HPFfiltCST);
                fluctCST = filtCST(tempSteadyEMGInd(1):tempSteadyEMGInd(2)-1);
                normFluctCST = fluctCST/max(abs(fluctCST));
                fineDelay = finddelay(normFluctForce,normFluctCST,maxDelay);
                fluctCST = filtCST(tempSteadyEMGInd(1)+fineDelay:tempSteadyEMGInd(2)+fineDelay-1);
                normFluctCST = fluctCST/max(abs(fluctCST));
                rho = corr(normFluctForce,normFluctCST');
                EMdelay = EMdelay - fineDelay;
                
                %save to excel     
                filename = "test.xlsx";
                pathname = "E:/Vibration";
                filename = fullfile(pathname,filename);
%                 EMdelayString = 'G';
%                 rhoString = 'J';
%                 fcHPFString = 'E';
%                 fcLPFString = 'F';
                EMdelayString = 'AU';
                rhoString = 'AX';
                fcHPFString = 'AS';
                fcLPFString = 'AT';
                writematrix(EMdelay,filename,'Sheet','Correlations','Range',EMdelayString+string(f+1));
                writematrix(rho,filename,'Sheet','Correlations','Range',rhoString+string(f+1));
                writematrix(fcHPF,filename,'Sheet','Correlations','Range',fcHPFString+string(f+1));
                writematrix(fcLPF,filename,'Sheet','Correlations','Range',fcLPFString+string(f+1));
                
                %%% populate excel for [fcHPF - fcLPF] force bandwidth %%%
                %read in EMdelays from excel
                EMdelay = EMdelaySamples(f);
                tempSteadyEMGInd = steadyEMGInd - EMdelay;
                
                %%% find optimal cutoff freqs %%%
                maxDelay = 0.2*fsamp; %max of 200 ms correction
                delayArray = zeros(6,9); %row 1,2=fcHPF 0.5Hz; row 3,4=fcHPF 0.75Hz; row 5,6=fcHPF 1Hz;
                                         %column=fcLPF 1, 1.5, 2, ..., 5 Hz                        
                i=1;
                for fcHPF = 0.5:0.25:1
                    [bHPF,aHPF] = butter(2,fcHPF/(fsamp/2),'high');
                    HPFfiltForce = filtfilt(bHPF,aHPF,force);
                    HPFfiltCST = filtfilt(bHPF,aHPF,smoothCST);
                    j=1;
                    for fcLPF = 1:0.5:5
                        [bLPF,aLPF] = butter(2,fcLPF/(fsamp/2));
                        [bLPFforce,aLPFforce] = butter(2,fcLPF/(fsamp/2));
                        filtForce = filtfilt(bLPFforce,aLPFforce,HPFfiltForce);
                        fluctForce = filtForce(steadyForceInd(1):steadyForceInd(2)-1);
                        normFluctForce = fluctForce/max(abs(fluctForce));
                        filtCST = filtfilt(bLPF,aLPF,HPFfiltCST);
                        fluctCST = filtCST(tempSteadyEMGInd(1):tempSteadyEMGInd(2)-1);
                        normFluctCST = fluctCST/max(abs(fluctCST));
                        fineDelay = finddelay(normFluctForce,normFluctCST,maxDelay);
                        fluctCST = filtCST(tempSteadyEMGInd(1)+fineDelay:tempSteadyEMGInd(2)+fineDelay-1);
                        normFluctCST = fluctCST/max(abs(fluctCST));
                        delayArray(1+2*(i-1),j) = corr(normFluctForce,normFluctCST');
                        delayArray(2+2*(i-1),j) = fineDelay;
                        j=j+1;
                    end
                    i=i+1;
                end
                rhoArray = [delayArray(1,:); delayArray(3,:); delayArray(5,:)];
                delayArray = [delayArray(2,:); delayArray(4,:); delayArray(6,:)];
                [rho,maxInd] = max(rhoArray,[],'all','linear');
                [maxInd1,maxInd2] = ind2sub([3,9],maxInd);
                fineDelay = delayArray(maxInd1,maxInd2);
                EMdelay = EMdelay - fineDelay;
                tempSteadyEMGInd = tempSteadyEMGInd + fineDelay;
                switch maxInd1
                    case 1
                        fcHPF = 0.5;
                    case 2
                        fcHPF = 0.75;
                    case 3
                        fcHPF = 1;
                end
                fcLPF = maxInd2*0.5 + 0.5;
                %save to excel     
                EMdelayString = 'AA';
                rhoString = 'AD';
                fcHPFString = 'Y';
                fcLPFString = 'Z';
                writematrix(EMdelay,filename,'Sheet','Correlations','Range',EMdelayString+string(f+1));
                writematrix(rho,filename,'Sheet','Correlations','Range',rhoString+string(f+1));
                writematrix(fcHPF,filename,'Sheet','Correlations','Range',fcHPFString+string(f+1));
                writematrix(fcLPF,filename,'Sheet','Correlations','Range',fcLPFString+string(f+1));
                
            elseif (run == 3)
                %read in EMdelays from excel
                EMdelay = EMdelaySamples(f);
                steadyEMGInd = steadyEMGInd - EMdelay;
            end
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% extract steady sections %%%
            %HPF smoothed CST at 0.5 Hz
            fcCST = 0.5;
            [b05, a05] = butter(2,fcCST/(fsamp/2),'high');
            steadyCST = filtfilt(b05,a05,smoothCST);
            steadyForce = force(steadyForceInd(1):steadyForceInd(2)-1);
            steadyCST = steadyCST(steadyEMGInd(1):steadyEMGInd(2)-1);
            steadyBinary = binary(:,steadyEMGInd(1):steadyEMGInd(2)-1);
            steadyActiveIndex = activeind(steadyEMGInd(1):steadyEMGInd(2)-1);
            steadyIPTs = IPTs(:,steadyEMGInd(1):steadyEMGInd(2)-1);

            %%% calculate PNRs over steady section %%%
            for mu=1:numMUs
                pulse=0; noise=0; p=0; n=0;
                for samp=2:numSteadyEMGSamps-1
                    if(steadyBinary(mu,samp)==1)
                        pulse = pulse + steadyIPTs(mu,samp)^2;
                        p=p+1;
                    elseif(steadyBinary(mu,samp-1)~=1 && steadyBinary(mu,samp+1)~=1)
                        noise = noise + steadyIPTs(mu,samp)^2;
                        n=n+1;
                    end
                end
                steadyPNRs(mu) = 10*log10((pulse/p)/(noise/n));
            end

            %%% calculate CVs over steady section %%%
            for win=1:numEpochs 
                winStartInd = stepSize*fsamp*(win-1) + 1; %1, 501, 1001, ..., 18001
                winEndInd = winStartInd + epochSize*fsamp - 1; %2000, 2500, 3000, ..., 20000
                tempForce = steadyForce(winStartInd:winEndInd);
                tempCST = steadyCST(winStartInd:winEndInd);
                CVforce(win) = 100*std(tempForce)/mean(tempForce);
                SDforce(win) = std(tempForce);
                SDCST(win) = std(tempCST);
                usableMUs(win) = numMUs;
                for mu=1:numMUs
                    tempBinary = steadyBinary(mu,winStartInd:winEndInd);
                    numISIs = sum(tempBinary)-1;
                    tempISIs = zeros(1,numISIs);
                    tempInd = find(tempBinary)+winStartInd-1;
                    for isi = 1:numISIs
                        tempISIs(isi) = (tempInd(isi+1)-tempInd(isi))/fsamp;
                    end
                    maxISIsamps = 333; %333 samples = 0.1665 s = ~6 Hz
                    if isempty(tempInd)
                        CVISI(mu,win) = NaN;
                        usableMUs(win) = usableMUs(win) - 1;
                    elseif (tempInd(1)-winStartInd > maxISIsamps || ...
                            winEndInd-tempInd(end) > maxISIsamps || ...
                            max(tempISIs) > maxISIsamps)
                        CVISI(mu,win) = NaN;
                        usableMUs(win) = usableMUs(win) - 1;
                    else
                        CVISI(mu,win) = 100*std(tempISIs)/mean(tempISIs);
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% populate structure                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            targetForce = append('t',string(t));
            VIBdata.subject(s).data.(conditionString).(targetForce).force = force;
            VIBdata.subject(s).data.(conditionString).(targetForce).emg = SIG;
            VIBdata.subject(s).data.(conditionString).(targetForce).cst = CST;
            VIBdata.subject(s).data.(conditionString).(targetForce).normcst = normCST;
            VIBdata.subject(s).data.(conditionString).(targetForce).binary = binary;
            VIBdata.subject(s).data.(conditionString).(targetForce).activeIndex = activeind;
            VIBdata.subject(s).data.(conditionString).(targetForce).fsamp = fsamp;
            VIBdata.subject(s).data.(conditionString).(targetForce).ISIs = ISIs;
            VIBdata.subject(s).data.(conditionString).(targetForce).MUPulses = MUPulses;
            VIBdata.subject(s).data.(conditionString).(targetForce).IPTs = IPTs;
            VIBdata.subject(s).data.(conditionString).(targetForce).PNRs = PNR;
            VIBdata.subject(s).data.(conditionString).(targetForce).steadyForceInd = steadyForceInd;
            VIBdata.subject(s).data.(conditionString).(targetForce).steadyEMGInd = steadyEMGInd;
            VIBdata.subject(s).data.(conditionString).(targetForce).EMdelay = EMdelay;
            VIBdata.subject(s).data.(conditionString).(targetForce).steadyForce = steadyForce;
            VIBdata.subject(s).data.(conditionString).(targetForce).steadyCST = steadyCST;
            VIBdata.subject(s).data.(conditionString).(targetForce).steadyBinary = steadyBinary;
            VIBdata.subject(s).data.(conditionString).(targetForce).steadyActiveIndex = steadyActiveIndex;
            VIBdata.subject(s).data.(conditionString).(targetForce).steadyIPTs = steadyIPTs;
            VIBdata.subject(s).data.(conditionString).(targetForce).steadyPNRs = steadyPNRs;
            VIBdata.subject(s).data.(conditionString).(targetForce).CVforce = CVforce;
            VIBdata.subject(s).data.(conditionString).(targetForce).SDforce = SDforce;
            VIBdata.subject(s).data.(conditionString).(targetForce).SDCST = SDCST;
            VIBdata.subject(s).data.(conditionString).(targetForce).CVISI = CVISI;
            VIBdata.subject(s).data.(conditionString).(targetForce).usableMUs = usableMUs;
            VIBdata.subject(s).data.(conditionString).(targetForce).fc = [fcHPF,fcLPF];
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%             %%% populate concatenated 5 and 20 percent %%%
%             if(t==20)
%                 % load variables into workspace
%                 filename = sprintf('%s-%s-5-20-Final.mat',VIBdata.subject(s).data.id,conditionString);
%                 pathname = 'E:/Vibration/Tracking/5 20/';
%                 filename = fullfile(pathname,filename);
%                 if isfile(filename)
%                     load(filename);     
%                     % calculate emg variables
%                     numMUs = numel(MUPulses);
%                     numSamples = floor(fsamp*SIGlength);
%                     binary = zeros(numMUs,numSamples);
%                     force = [VIBdata.subject(s).data.(conditionString).t5.force; ...
%                         VIBdata.subject(s).data.(conditionString).t20.force];
%                     
%                     % populate binary
%                     for mu = 1:numMUs
%                         for spikeIndex = 1:numel(MUPulses{mu})
%                             tempSpikeTime = MUPulses{mu}(spikeIndex);
%                             binary(mu,tempSpikeTime) = 1;
%                         end
%                     end
% 
%                     % calculate ISIs
%                     ISIs = cell(1,size(MUPulses,2));
%                     for mu = 1:numMUs
%                         tempISI = zeros(1,numel(MUPulses{mu})-1);
%                         for spikeIndex = 2:numel(MUPulses{mu})
%                             tempISI(1,spikeIndex-1) = MUPulses{mu}(spikeIndex)-MUPulses{mu}(spikeIndex-1);
%                         end
%                         ISIs{mu} = tempISI;
%                     end
% 
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %%% calculate variables over steady 10 seconds              %%%
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     steadyEMGInd = cell(2,1);
%                     CVISI = cell(2,1);
%                     steadyPNRs = cell(2,1);
%                     steadyBinary = cell(2,1);
%                     steadyIPTs = cell(2,1);
% 
%                     %cell array variables, row 1 = 5%, row 2 = 20%
%                     for i=1:2
%                         targetForce = sprintf('t%d',5+(i-1)*15);
%                         switch i
%                             case 1
%                                 steadyEMGInd{i,1} = VIBdata.subject(s).data.(conditionString).(targetForce).steadyEMGInd;
%                             case 2
%                                 steadyEMGInd{i,1} = VIBdata.subject(s).data.(conditionString).(targetForce).steadyEMGInd + ...
%                                                         length(VIBdata.subject(s).data.(conditionString).t5.force);
%                         end
%                         numSteadyEMGSamps = 10*fsamp;
%                         CVISI{i,1} = zeros(numMUs,numEpochs);
%                         steadyPNRs{i,1} = zeros(1,numMUs);
%                         steadyBinary{i,1} = binary(:,steadyEMGInd{i,1}(1):steadyEMGInd{i,1}(2)-1);
%                         steadyIPTs{i,1} = IPTs(:,steadyEMGInd{i,1}(1):steadyEMGInd{i,1}(2)-1);
% 
%                         %%% calculate PNRs over steady section %%%
%                         for mu=1:numMUs
%                             pulse=0; noise=0; p=0; n=0;
%                             for samp=2:numSteadyEMGSamps-1
%                                 if(steadyBinary{i,1}(mu,samp)==1)
%                                     pulse = pulse + steadyIPTs{i,1}(mu,samp)^2;
%                                     p=p+1;
%                                 elseif(steadyBinary{i,1}(mu,samp-1)~=1 && steadyBinary{i,1}(mu,samp+1)~=1)
%                                     noise = noise + steadyIPTs{i,1}(mu,samp)^2;
%                                     n=n+1;
%                                 end
%                             end
%                             steadyPNRs{i,1}(mu) = 10*log10((pulse/p)/(noise/n));
%                         end
% 
%                         %%% calculate CVISIs over steady section - 1s windows %%%
%                         for win=1:numEpochs % steady 10s, 1s sliding window, 0.25s step size
%                             winStartInd = 0.25*fsamp*(win-1) + 1; %1, 501, 1001, ..., 18001
%                             winEndInd = winStartInd + fsamp - 1; %2000, 2500, 3000, ..., 20000
%                             for mu=1:numMUs
%                                 tempBinary = steadyBinary{i,1}(mu,winStartInd:winEndInd);
%                                 numISIs = sum(tempBinary)-1;
%                                 tempISIs = zeros(1,numISIs);
%                                 tempInd = find(tempBinary)+winStartInd-1;
%                                 for isi = 1:numISIs
%                                     tempISIs(isi) = (tempInd(isi+1)-tempInd(isi))/fsamp;
%                                 end
%                                 maxISIsamps = 333; %333 samples = 0.1665 s = ~6 Hz
%                                 if isempty(tempInd)
%                                     CVISI{i,1}(mu,win) = NaN;
%                                 elseif (tempInd(1)-winStartInd > maxISIsamps || ...
%                                         winEndInd-tempInd(end) > maxISIsamps || ...
%                                         max(tempISIs) > maxISIsamps)
%                                     CVISI{i,1}(mu,win) = NaN;
%                                 else
%                                     CVISI{i,1}(mu,win) = 100*std(tempISIs)/mean(tempISIs);
%                                 end
%                             end
%                         end
%                     end
% 
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %%% populate structure                                      %%%
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     VIBdata.subject(s).data.(conditionString).t520.force = force;
%                     VIBdata.subject(s).data.(conditionString).t520.emg = SIG;
%                     VIBdata.subject(s).data.(conditionString).t520.binary = binary;
%                     VIBdata.subject(s).data.(conditionString).t520.fsamp = fsamp;
%                     VIBdata.subject(s).data.(conditionString).t520.ISIs = ISIs;
%                     VIBdata.subject(s).data.(conditionString).t520.MUPulses = MUPulses;
%                     VIBdata.subject(s).data.(conditionString).t520.IPTs = IPTs;
%                     VIBdata.subject(s).data.(conditionString).t520.PNRs = PNR;
%                     VIBdata.subject(s).data.(conditionString).t520.steadyEMGInd = steadyEMGInd;
%                     VIBdata.subject(s).data.(conditionString).t520.steadyBinary = steadyBinary;
%                     VIBdata.subject(s).data.(conditionString).t520.steadyIPTs = steadyIPTs;
%                     VIBdata.subject(s).data.(conditionString).t520.steadyPNRs = steadyPNRs;
%                     VIBdata.subject(s).data.(conditionString).t520.CVISI = CVISI;
%                 end
%             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%             %%% populate concatenated vision, sound, vibration on last iteration %%%
%             if(c==3 && t==20)
%                 for j=1:2 %loop through 5 and 20 percent
%                     targetForce = sprintf('t%d',5+(j-1)*15);
%                     % load variables into workspace from concatenated file
%                     filename = sprintf('%s_Post_edited.mat',VIBdata.subject(s).data.id);
%                     pathname = 'E:/Synaptic Noise/FDI Edited/5 15 25/';
%                     filename = fullfile(pathname,filename);
%                     load(filename);     
% 
%                     numSamples = floor(fsamp*SIGlength);
%                     CST = zeros(1,numSamples);
%                     binary = zeros(size(MUPulses,2),numSamples);
%                     start = MUPulses{1}(1);
%                     endd = MUPulses{1}(end);
% 
%                     % populate CST
%                     for mu = 1:numel(MUPulses)
%                         for spikeIndex = 1:numel(MUPulses{mu})
%                             tempSpikeTime = MUPulses{mu}(spikeIndex);
%                             binary(mu,tempSpikeTime) = 1;
%                             CST(tempSpikeTime) = CST(tempSpikeTime) + 1;
%                             if MUPulses{mu}(spikeIndex) < start
%                                 start = MUPulses{mu}(spikeIndex);
%                             end
%                             if MUPulses{mu}(spikeIndex) > endd
%                                 endd = MUPulses{mu}(spikeIndex);
%                             end
%                         end
%                     end
% 
%                     % calculate active index and normalize CST
%                     win = 0.25*fsamp; %250 ms window
%                     activeind = [];
%                     for r = 1:size(binary,1)
%                         start = 1;
%                         endd = size(binary,2);
%                         num = floor((endd-start)/win)-1;
%                         for i = 0:num-1 % slide 200 ms window to calculate active index
%                             ind = start+(win*i);
%                             temp = binary(r,ind:ind+win);
%                             if sum(temp) > 0
%                                 activeind(r,ind:ind+win) = 1;
%                             elseif sum(temp) == 0
%                                 activeind(r,ind:ind+win) = 0;
%                             end
%                         end
%                     end
% 
%                     buff = length(binary)-length(activeind(1,:));
%                     zer = zeros(size(activeind,1),buff);
%                     newAI= horzcat(activeind,zer);
%                     numMU = sum(newAI);
%                     activeind = numMU;
%                     normCST = CST./activeind;
%                     normCST(isinf(normCST)|isnan(normCST)) = 0;
% 
%                     % calculate ISIs
%                     ISIs = cell(1,size(MUPulses,2));
%                     for mu = 1:numel(MUPulses)
%                         tempISI = zeros(1,numel(MUPulses{mu})-1);
%                         for spikeIndex = 2:numel(MUPulses{mu})
%                             tempISI(1,spikeIndex-1) = MUPulses{mu}(spikeIndex)-MUPulses{mu}(spikeIndex-1);
%                         end
%                         ISIs{mu} = tempISI;
%                     end
%                 
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %%% calculate variables over steady 10 seconds              %%%
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     numMUs = numel(MUPulses);
% 
%                     %%% smooth and filter entire signals %%%
%                     %smooth CST
%                     hann_window = hann(floor(0.4*fsamp));
%                     smoothCST = conv(normCST,hann_window,'same');
% 
%                     %cell array vars, row1 = vision, row2 = sound, row3 = vibration
%                     for i=1:3
%                         steadyEMGInd{i,1} = VIBdata.subject(s).data.condition(i).(targetForce).steadyEMGInd;
%                         numSteadyEMGSamps = 10*fsamp;
%                         CVISI{i,1} = zeros(numMUs,10);
%                         SDCST{i,1} = zeros(1,10);
%                         steadyPNRs{i,1} = zeros(1,numMUs);
% 
%                         %%% extract steady sections %%%
%                         %HPF smoothed CST at 0.5 Hz
%                         fcCST = 0.5;
%                         steadyCST{i,1} = filtfilt(b05,a05,smoothCST);
%                         steadyCST{i,1} = steadyCST{i,1}(steadyEMGInd{i,1}(1):steadyEMGInd{i,1}(2)-1);
%                         steadyActiveIndex{i,1} = activeind(steadyEMGInd{i,1}(1):steadyEMGInd{i,1}(2)-1);
%                         steadyBinary{i,1} = binary(:,steadyEMGInd{i,1}(1):steadyEMGInd{i,1}(2)-1);
%                         steadyIPTs{i,1} = IPTs(:,steadyEMGInd{i,1}(1):steadyEMGInd{i,1}(2)-1);
% 
%                         %%% calculate PNRs over steady section %%%
%                         for mu=1:numMUs
%                             pulse=0; noise=0; p=0; n=0;
%                             for samp=2:numSteadyEMGSamps-1
%                                 if(steadyBinary{i,1}(mu,samp)==1)
%                                     pulse = pulse + steadyIPTs{i,1}(mu,samp)^2;
%                                     p=p+1;
%                                 elseif(steadyBinary{i,1}(mu,samp-1)~=1 && steadyBinary{i,1}(mu,samp+1)~=1)
%                                     noise = noise + steadyIPTs{i,1}(mu,samp)^2;
%                                     n=n+1;
%                                 end
%                             end
%                             steadyPNRs{i,1}(mu) = 10*log10((pulse/p)/(noise/n));
%                         end
% 
%                         %%% calculate CVs over steady section - 1s windows %%%
%                         for sec=1:10
%                             emgWindowStart = 1 + fsamp*(sec-1);
%                             tempCST = steadyCST{i,1}(emgWindowStart:emgWindowStart+fsamp-1);
%                             SDCST{i,1}(sec) = std(tempCST);
%                             for mu=1:numMUs
%                                 tempBinary = steadyBinary{i,1}(mu,emgWindowStart:emgWindowStart+fsamp-1);
%                                 %need at least 3 spikes (2 ISIs) to calculate CV
%                                 if(sum(tempBinary)<3)
%                                     CVISI{i,1}(mu,sec) = NaN;
%                                 else
%                                     numISIs = sum(tempBinary)-1;
%                                     tempISIs = zeros(1,numISIs);
%                                     tempInd = find(tempBinary);
%                                     for isi=1:numISIs
%                                         tempISIs(isi) = fsamp*(tempInd(isi+1)-tempInd(isi));
%                                     end
%                                     CVISI{i,1}(mu,sec) = 100*std(tempISIs)/mean(tempISIs);
%                                 end
%                             end
%                         end
%                     end
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %%% populate structure                                      %%%
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     VIBdata.subject(s).data.concatenated.(targetForce).emg = SIG;
%                     VIBdata.subject(s).data.concatenated.(targetForce).cst = CST;
%                     VIBdata.subject(s).data.concatenated.(targetForce).normcst = normCST;
%                     VIBdata.subject(s).data.concatenated.(targetForce).binary = binary;
%                     VIBdata.subject(s).data.concatenated.(targetForce).activeIndex = activeind;
%                     VIBdata.subject(s).data.concatenated.(targetForce).fsamp = fsamp;
%                     VIBdata.subject(s).data.concatenated.(targetForce).ISIs = ISIs;
%                     VIBdata.subject(s).data.concatenated.(targetForce).MUPulses = MUPulses;
%                     VIBdata.subject(s).data.concatenated.(targetForce).IPTs = IPTs;
%                     VIBdata.subject(s).data.concatenated.(targetForce).PNRs = PNR;
%                     VIBdata.subject(s).data.concatenated.(targetForce).steadyEMGInd = steadyEMGInd;
%                     VIBdata.subject(s).data.concatenated.(targetForce).steadyCST = steadyCST;
%                     VIBdata.subject(s).data.concatenated.(targetForce).steadyBinary = steadyBinary;
%                     VIBdata.subject(s).data.concatenated.(targetForce).steadyActiveIndex = steadyActiveIndex;
%                     VIBdata.subject(s).data.concatenated.(targetForce).steadyIPTs = steadyIPTs;
%                     VIBdata.subject(s).data.concatenated.(targetForce).steadyPNRs = steadyPNRs;
%                     VIBdata.subject(s).data.concatenated.(targetForce).SDCST = SDCST;
%                     VIBdata.subject(s).data.concatenated.(targetForce).CVISI = CVISI;
%                 end
%             end
            
            % increase force file index
            f=f+1;
        end
    end
    clearvars -except VIBdata fsamp CVStartIndex f EMdelaySamples TrialNum run fc forceLPF ...
        numEpochs epochSize stepSize;
    fclose('all');
end

% save('VIBdata.mat','VIBdata','-v7.3');
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

























