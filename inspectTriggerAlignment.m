clearvars -except VIBdata;
clc;

TrialNum = double(readmatrix('E:/Vibration/Participant_data.xlsx','Sheet','General','Range','D2:D169'));
f=1;

for s = 1:length(VIBdata.subject)
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
        for t=[5,20]
            %import vars
            targetForce = append('t',string(t));
            subjectID = VIBdata.subject(s).data.id;
            force = VIBdata.subject(s).data.(conditionString).(targetForce).force;
            emg = VIBdata.subject(s).data.(conditionString).(targetForce).emg;
            normCST = VIBdata.subject(s).data.(conditionString).(targetForce).normcst;
            fsamp = VIBdata.subject(s).data.(conditionString).(targetForce).fsamp;
            MUPulses = VIBdata.subject(s).data.(conditionString).(targetForce).MUPulses;
            steadyForceInd = VIBdata.subject(s).data.(conditionString).(targetForce).steadyForceInd;
            steadyEMGInd = VIBdata.subject(s).data.(conditionString).(targetForce).steadyEMGInd;
            EMdelay = VIBdata.subject(s).data.(conditionString).(targetForce).EMdelay;
            steadyBinary = VIBdata.subject(s).data.(conditionString).(targetForce).steadyBinary;
            numMUs = size(steadyBinary,1);
            
            %%% smooth and filter entire signals %%%
            %smooth CST
            hann_window = hann(floor(0.4*fsamp));
            smoothCST = conv(normCST,hann_window,'same');
            
            %filter signals
            fcHPF = 1;
            fcLPF = 2;
            forceLPF = 2;
            [bHPF,aHPF] = butter(2,fcHPF/(fsamp/2),'high');
            [bLPF,aLPF] = butter(2,fcLPF/(fsamp/2));
            [bLPFforce,aLPFforce] = butter(2,forceLPF/(fsamp/2));
            filtCST = filtfilt(bHPF,aHPF,smoothCST);
            filtCST = filtfilt(bLPF,aLPF,filtCST);
            filtForce = filtfilt(bHPF,aHPF,force);
            filtForce = filtfilt(bLPFforce,aLPFforce,filtForce);
            fluctForce = filtForce(steadyForceInd(1):steadyForceInd(2)-1);
            fluctCST = filtCST(steadyEMGInd(1):steadyEMGInd(2)-1);
            fluctCST2 = filtCST(steadyForceInd(1):steadyForceInd(2)-1);
            iEMG = double(emg{8,4})';
            iEMG = bandpass(iEMG,[20 500],fsamp);
            iEMG = abs(iEMG);
            [bEMG,aEMG] = butter(2,2/(fsamp/2));
            iEMG = filtfilt(bEMG,aEMG,iEMG);
            
            %%% plot %%%
            steadyForceOffset = double(steadyForceInd(1))/fsamp;
            time = 0:1/fsamp:(length(force)-1)/fsamp;
            time10 = 0:1/fsamp:10-1/fsamp; 
            EMdelayMS = round(1000*EMdelay/fsamp); 
            normFluctForce = fluctForce/max(abs(fluctForce)); 
            normFluctCST = fluctCST/max(abs(fluctCST)); 
            normFluctCST2 = fluctCST2/max(abs(fluctCST2));
            rho = round(corr(normFluctForce,normFluctCST'),2);
            rho2 = round(corr(normFluctForce,normFluctCST2'),2);
            fig = figure(1);
            
            %force and EMG envelopes
            subplot(2,1,1);
            hold on;
            plot(time,force/max(abs(force)),...
                time10+steadyForceOffset,force(steadyForceInd(1):steadyForceInd(2)-1)/max(abs(force)),...
                time+EMdelay/fsamp,iEMG/max(abs(iEMG)),...
                time,iEMG/max(abs(iEMG)),...
                'LineWidth',2);
            steadyLabels = [sprintf("x = %0.2f",steadyForceOffset) sprintf("x = %0.2f",steadyForceOffset+10)];
            xline([steadyForceOffset steadyForceOffset+10],'-',{steadyLabels(1),steadyLabels(2)},...
                'LabelVerticalAlignment','bottom');
            legend('force','steady10','manual align','trigger align'); 
            xlabel("Time (s)");
            title("Force and EMG Envelopes"); 
            hold off;

            %force/cst fluctuations
            subplot(2,1,2); 
            hold on;
            plot(time10+steadyForceOffset,normFluctForce,'LineWidth',2,'Color',"#0072BD");
            plot(time10+steadyForceOffset,normFluctCST,'LineWidth',2,'Color',"#EDB120"); 
            plot(time10+steadyForceOffset,normFluctCST2,'LineWidth',2,'Color',"#7E2F8E");
            xlim('tight');
            legendString = sprintf("CST [%0.2f - %0.1f Hz] (r = %0.2f)",fcHPF,fcLPF,rho);
            legendString2 = sprintf("CST [%0.2f - %0.1f Hz] (r = %0.2f)",fcHPF,fcLPF,rho2);
            legend('force',legendString,legendString2);
            xlabel("Time (s)");
            title("Force and CST Fluctuations (steady 10s)");
            hold off;

            titleString = sprintf("%s %s %d%% %d MUs (%d ms difference)",VIBdata.subject(s).data.id,conditionString,t,numel(MUPulses),EMdelayMS);
            sgtitle(titleString);
            datacursormode on;    
            waitforbuttonpress;    
            clf(fig);
            clearvars -except VIBdata f s c TrialNum conditionString
            f=f+1;
        end
    end
end