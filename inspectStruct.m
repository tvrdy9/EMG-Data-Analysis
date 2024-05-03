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
            CST = VIBdata.subject(s).data.(conditionString).(targetForce).cst;
            normCST = VIBdata.subject(s).data.(conditionString).(targetForce).normcst;
            binary = VIBdata.subject(s).data.(conditionString).(targetForce).binary;
            activeIndex = VIBdata.subject(s).data.(conditionString).(targetForce).activeIndex;
            fsamp = VIBdata.subject(s).data.(conditionString).(targetForce).fsamp;
            ISIs = VIBdata.subject(s).data.(conditionString).(targetForce).ISIs;
            MUPulses = VIBdata.subject(s).data.(conditionString).(targetForce).MUPulses;
            IPTs = VIBdata.subject(s).data.(conditionString).(targetForce).IPTs;
            PNRs = VIBdata.subject(s).data.(conditionString).(targetForce).PNRs;
            steadyForceInd = VIBdata.subject(s).data.(conditionString).(targetForce).steadyForceInd;
            steadyEMGInd = VIBdata.subject(s).data.(conditionString).(targetForce).steadyEMGInd;
            EMdelay = VIBdata.subject(s).data.(conditionString).(targetForce).EMdelay;
            steadyForce = VIBdata.subject(s).data.(conditionString).(targetForce).steadyForce;
            filtsmoothnormCST = VIBdata.subject(s).data.(conditionString).(targetForce).steadyCST;
            steadyBinary = VIBdata.subject(s).data.(conditionString).(targetForce).steadyBinary;
            steadyActiveIndex = VIBdata.subject(s).data.(conditionString).(targetForce).steadyActiveIndex;
            steadyIPTs = VIBdata.subject(s).data.(conditionString).(targetForce).steadyIPTs;
            steadyPNRs = VIBdata.subject(s).data.(conditionString).(targetForce).steadyPNRs;
            CVforce = VIBdata.subject(s).data.(conditionString).(targetForce).CVforce;
            SDCST = VIBdata.subject(s).data.(conditionString).(targetForce).SDCST;
            CVISI = VIBdata.subject(s).data.(conditionString).(targetForce).CVISI;
            usableMUs = VIBdata.subject(s).data.(conditionString).(targetForce).usableMUs;
            fc = VIBdata.subject(s).data.(conditionString).(targetForce).fc;
            numMUs = size(steadyBinary,1);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %%% inspect steady 10s %%%
            %%% smooth and filter entire signals %%%
            %smooth CST
            hann_window = hann(floor(0.4*fsamp));
            smoothCST = conv(normCST,hann_window,'same');
            
            %filter signals
            fcHPF = fc(1);
            fcLPF = fc(2);
            forceLPF = 5;
            [bHPF,aHPF] = butter(2,fcHPF/(fsamp/2),'high');
            [bLPF,aLPF] = butter(2,fcLPF/(fsamp/2));
            [bLPFforce,aLPFforce] = butter(2,forceLPF/(fsamp/2));
            filtCST = filtfilt(bHPF,aHPF,smoothCST);
            filtCST = filtfilt(bLPF,aLPF,filtCST);
            filtForce = filtfilt(bHPF,aHPF,force);
            filtForce = filtfilt(bLPFforce,aLPFforce,filtForce);
            fluctForce = filtForce(steadyForceInd(1):steadyForceInd(2)-1);
            fluctCST = filtCST(steadyEMGInd(1):steadyEMGInd(2)-1);
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
            rho = round(corr(normFluctForce,normFluctCST'),2);
            fig = figure(1);
            %save p-values
%             cellString = sprintf('L%d',f+1);
%             filename = "Participant_data.xlsx";
%             pathname = "D:/Vibration";
%             filename = fullfile(pathname,filename);
%             writematrix(rho,filename,'Sheet','General','Range',cellString);
            
            %force and EMG envelopes
            subplot(3,1,1);
            hold on;
            plot(time,force/max(abs(force)),...
                time10+steadyForceOffset,force(steadyForceInd(1):steadyForceInd(2)-1)/max(abs(force)),...
                time+EMdelay/fsamp,iEMG/max(abs(iEMG)),...
                'LineWidth',2);
            steadyLabels = [sprintf("x = %0.2f",steadyForceOffset) sprintf("x = %0.2f",steadyForceOffset+10)];
            xline([steadyForceOffset steadyForceOffset+10],'-',{steadyLabels(1),steadyLabels(2)},...
                'LabelVerticalAlignment','bottom');
            legend('force','steady10','emg'); 
            xlabel("Time (s)");
            title("Force and EMG Envelopes"); 
            hold off;

            %force/cst fluctuations
            subplot(3,1,2); 
            hold on;
            plot(time10+steadyForceOffset,normFluctForce,'LineWidth',2,'Color',"#0072BD");
            plot(time10+steadyForceOffset,normFluctCST,'LineWidth',2,'Color',"#EDB120"); 
            xlim('tight');
            legendString = sprintf("CST [%0.2f - %0.1f Hz] (r=%0.2f)",fcHPF,fcLPF,rho);
            legend('force',legendString);
            xlabel("Time (s)");
            title("Force and CST Fluctuations");
            hold off;

            %discharge characteristics
            subplot(3,1,3);
            hold on;
            labelArray = strings(numMUs,1);
            badPNRs = zeros(numMUs,1);
            %create raster plot
            for mu=1:numMUs
                tempSpikeIndeces = find(steadyBinary(mu,:));
                tempSpikeTimes = tempSpikeIndeces/fsamp + steadyForceOffset;
                if isnan(PNRs(mu))
                    PNRs(mu) = 0;
                elseif isnan(steadyPNRs(mu))
                    steadyPNRs(mu) = 0;
                end
                for spike = tempSpikeTimes
                    line([spike spike],[mu-0.5 mu+0.5],'Color','blue');
                    labelArray(mu) = sprintf("MU%d PNR = %0.2f dB (%0.2f dB)",mu,PNRs(mu),steadyPNRs(mu));
                    if (PNRs(mu) < 30 && steadyPNRs(mu) < 30)
                        badPNRs(mu) = 1;
                    end
                end
            end
            xlim('tight');
            ylim([0.5, numMUs+0.5]);
            xline(steadyForceOffset,'k','LineWidth',2);
            for i=1:40
                xline(steadyForceOffset+i*0.25,'k','LineWidth',2);
            end
            %plot red spikes/epochs
            for i=1:37
                for mu=1:numMUs
                    if isnan(CVISI(mu,i))
                        tempSpikeIndeces = find(steadyBinary(mu,1+(i-1)*0.25*fsamp:fsamp+(i-1)*0.25*fsamp));
                        tempSpikeTimes = tempSpikeIndeces/fsamp + steadyForceOffset + (i-1)*0.25 - 1/fsamp;
                        for spike = tempSpikeTimes
                            line([spike spike],[mu-0.5 mu+0.5],'Color','red');
                        end
                        if usableMUs(i) < 5
                            xline(steadyForceOffset+(i-1)*0.25,'r','LineWidth',2);
                            xline(steadyForceOffset+(i-1)*0.25+1,'r','LineWidth',2);
                        end
                    end
                end
            end
            x = steadyForceOffset+0.025:0.25:steadyForceOffset+9.025;
            y = usableMUs + " MUs";
            for i=1:37
                if usableMUs(i) < 5
                    y{i} = ['\color{red}' y{i}];
                end
            end
            text(x,ones(1,37),y,'FontSize',8);
            xticks(steadyForceOffset:1:steadyForceOffset+10);
            xticklabels(0:1:10);
            xlabel("Epochs (s)");
            yticks(1:numMUs);
            labelArrayRed = cell(size(labelArray));
            for i=1:numMUs
                if badPNRs(i)
                    labelArrayRed{i} = ['\color{red}' labelArray{i}];
                else
                    labelArrayRed{i} = ['\color{black}' labelArray{i}];
                end
            end
            set(gca,'YTickLabel',labelArrayRed);
            title("Discharge Characteristics");
            titleString = sprintf("%s %s %d%% %d MUs (%d ms delay)",VIBdata.subject(s).data.id,conditionString,t,numel(MUPulses),EMdelayMS);
            sgtitle(titleString);
            datacursormode on;     
            hold off;
            waitforbuttonpress;    
            clf(fig);
            clearvars -except VIBdata f s c TrialNum conditionString
            f=f+1;
        end
    end
end