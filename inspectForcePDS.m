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
            numMUs = size(steadyBinary,1);

            L = length(force);
            PDS = fft(force); P2 = abs(PDS/L); P1 = P2(1:L/2+1); P1(2:end-1) = 2*P1(2:end-1);
            freq = fsamp*(0:(L/2))/L;
            plot(freq,P1); xlim([0 10]); ylim([0 0.1]);
            xlabel("Freq (Hz)"); ylabel("Power (au)"); title("Force PDS");
            hold on;
            f=f+1;
        end
    end
end
