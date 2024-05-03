clearvars -except VIBdata;
clc;

TrialNum = double(readmatrix('E:/Vibration/Participant_data.xlsx','Sheet','General','Range','D2:D169'));
f=1;
epochSize = 1;

% %0.75 - 5 Hz alignment
% EMdelaySamples = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','G2:G169'));
% fc = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','E2:F169'));
% sheetString = '0.75-5Hz';
% steadyCVISIavgString = 'H';
% steadySDForceString = 'I';
% epochCVForceSDCSTrhoString = 'M';
% epochCVForceCVISIavgString = 'S';
% epochSDForceSDCSTrhoString = 'P';
% epochSDForceCVISIavgString = 'V';
            
% %fcHPF - fcLPF Hz alignment
% EMdelaySamples = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','AA2:AA169'));
% fc = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','Y2:Z169'));
% sheetString = 'fcHPF-fcLPF';
% steadyCVISIavgString = 'AB';
% steadySDForceString = 'AC';
% epochCVForceSDCSTrhoString = 'AG';
% epochCVForceCVISIavgString = 'AM';
% epochSDForceSDCSTrhoString = 'AJ';
% epochSDForceCVISIavgString = 'AP';

%1.5 - 5 Hz alignment
EMdelaySamples = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','AU2:AU169'));
fc = double(readmatrix('E:/Vibration/test.xlsx','Sheet','Correlations','Range','AS2:AT169'));
sheetString = '1.5-5Hz';
steadyCVISIavgString = 'AV';
steadySDForceString = 'AW';
epochCVForceSDCSTrhoString = 'BA';
epochCVForceCVISIavgString = 'BG';
epochSDForceSDCSTrhoString = 'BD';
epochSDForceCVISIavgString = 'BJ';

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
            SDforce = VIBdata.subject(s).data.(conditionString).(targetForce).SDforce;
            SDCST = VIBdata.subject(s).data.(conditionString).(targetForce).SDCST;
            CVISI = VIBdata.subject(s).data.(conditionString).(targetForce).CVISI;
            usableMUs = VIBdata.subject(s).data.(conditionString).(targetForce).usableMUs;
            numMUs = size(steadyBinary,1);

            %%% Calculate correlations %%%
            steadySDForce = std(force(steadyForceInd(1):steadyForceInd(2)));
            CVForceCVISIrho = zeros(numMUs,1);
            SDForceCVISIrho = zeros(numMUs,1);
            epochCVForceSDCSTrho = corr(CVforce',SDCST','rows','complete');
            epochSDForceSDCSTrho = corr(SDforce',SDCST','rows','complete');
            steadyCVISIs = zeros(numMUs,1);
            for mu=1:numMUs
                tempBinary = steadyBinary(mu,:);
                numISIs = sum(tempBinary)-1;
                tempISIs = zeros(1,numISIs);
                tempInd = find(tempBinary);
                for isi = 1:numISIs
                    tempISIs(isi) = NaN;
                    if tempInd(isi+1)-tempInd(isi) <= 333
                        tempISIs(isi) = (tempInd(isi+1)-tempInd(isi))/fsamp;
                    end
                end
                steadyCVISIs(mu) = nanmean(tempISIs);
                if nnz(~isnan(CVISI(mu,:))) >= 5
                    CVForceCVISIrho(mu) = corr(CVforce',CVISI(mu,:)','rows','complete');
                    SDForceCVISIrho(mu) = corr(SDforce',CVISI(mu,:)','rows','complete');
                else
                    CVForceCVISIrho(mu) = NaN;
                    SDForceCVISIrho(mu) = NaN;
                end
            end
            steadyCVISIavg = mean(steadyCVISIs);
            epochCVForceCVISIavg = mean(CVForceCVISIrho,'omitnan');
            epochSDForceCVISIavg = mean(SDForceCVISIrho,'omitnan');
            
            switch epochSize
                case 1
                    CVforceString = sprintf('E%d:AO%d',f+2,f+2);
                    CVISIavgString = sprintf('AP%d:BZ%d',f+2,f+2);
                    SDCSTstring = sprintf('CA%d:DK%d',f+2,f+2);
                    SDforceString = sprintf('DL%d:EV%d',f+2,f+2);
                case 2
                    CVforceString = sprintf('E%d:U%d',f+2,f+2);
                    CVISIavgString = sprintf('V%d:AL%d',f+2,f+2);
                    SDCSTstring = sprintf('AM%d:BC%d',f+2,f+2);
                    SDforceString = sprintf('BD%d:BT%d',f+2,f+2);
                case 4
                    CVforceString = sprintf('E%d:K%d',f+2,f+2);
                    CVISIavgString = sprintf('L%d:R%d',f+2,f+2);
                    SDCSTstring = sprintf('S%d:Y%d',f+2,f+2);
                    SDforceString = sprintf('Z%d:AF%d',f+2,f+2);
            end
            
            % save vars to excel
            filename = "test.xlsx";
            pathname = "E:/Vibration";
            filename = fullfile(pathname,filename);
            writematrix(epochCVForceSDCSTrho,filename,'Sheet','Correlations','Range',epochCVForceSDCSTrhoString+string(f+1));
            writematrix(epochCVForceCVISIavg,filename,'Sheet','Correlations','Range',epochCVForceCVISIavgString+string(f+1));
            writematrix(epochSDForceCVISIavg,filename,'Sheet','Correlations','Range',epochSDForceCVISIavgString+string(f+1));
            writematrix(epochSDForceSDCSTrho,filename,'Sheet','Correlations','Range',epochSDForceSDCSTrhoString+string(f+1));            
            writematrix(steadyCVISIavg,filename,'Sheet','Correlations','Range',steadyCVISIavgString+string(f+1));
            writematrix(steadySDForce,filename,'Sheet','Correlations','Range',steadySDForceString+string(f+1));
            writematrix(CVforce,filename,'Sheet',sheetString,'Range',CVforceString);
            writematrix(mean(CVISI,'omitnan'),filename,'Sheet',sheetString,'Range',CVISIavgString);
            writematrix(SDCST,filename,'Sheet',sheetString,'Range',SDCSTstring);
            writematrix(SDforce,filename,'Sheet',sheetString,'Range',SDforceString);
            f=f+1;
        end
    end
end

        