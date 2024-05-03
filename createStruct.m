%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code builds a structure to hold all subject data with the                  %%%
%%% following architecture.                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                                 %%%
%%% VIBdata                           %structure containing all data                %%%
%%%     subject                       %structure array of subjects                  %%%
%%%         data                      %subject data                                 %%%
%%%             id                    %subject id code                              %%%
%%%             condition             %structure array of trial conditions          %%%
%%%                                   %     condition(1)=vision                     %%%
%%%                                   %     condition(2)=sound                      %%%
%%%                                   %     condition(3)=vibration                  %%%
%%%                 t5                %target force 5% MVC                          %%%
%%%                     force         %force signal trimmed after trigger, 20Hz LPF %%%
%%%                     emg           %emg signal trimmed around first ramp         %%%
%%%                     cst           %unedited cumulative spike train              %%%
%%%                     normcst       %normalized to 250ms activeIndex              %%%
%%%                     binary        %binary MU spike trains                       %%%
%%%                     activeIndex   %250ms window, 1=active, 0=inactive           %%%
%%%                     fsamp         %sampling freq of EMG                         %%%
%%%                     ISIs          %inter-spike interval (current-prev)          %%%
%%%                     MUPulses      %sample numbers of MU spikes                  %%%
%%%                     IPTs          %estimated spike trains                       %%%
%%%                     PNRs          %pulse-to-noise ratios for each MU            %%%
%%%                     steadyForceInd   %first and last sample #s                  %%%
%%%                     steadyEMGInd     %first and last sample #s                  %%%
%%%                     EMdelay       %electro-mechanical delay applied to EMG vars %%%
%%%                     steadyForce   %steady 10s force data                        %%%
%%%                     steadyCST     %steady 10s of normalized+smoothed CST        %%%
%%%                     steadyBinary  %steady 10s of binary MU spike trains         %%%
%%%                     steadyActiveIndex  %steady 10s of activeIndex               %%%
%%%                     steadyIPTs    %steady 10s of IPTs                           %%%
%%%                     steadyPNRs    %PNRs calculated over steady 10s              %%%
%%%                     CVforce       %CV for force from 1s epochs of steady 10s    %%%
%%%                     SDCST         %SD of CST from 1s epochs of steady 10s       %%%
%%%                     CVISI         %CV for ISI from 1s epochs of steady 10s      %%%
%%%                                                                                 %%%
%%%                 t20                                                             %%%
%%%                     ...                                                         %%%
%%%                                                                                 %%%
%%%                 t520              %both target forces concatenated              %%%
%%%                     emg           %emg signal trimmed after trigger             %%%
%%%                     cst           %unedited cumulative spike train              %%%
%%%                     normcst       %normalized to 250ms activeIndex              %%%
%%%                     binary        %binary MU spike trains                       %%%
%%%                     activeIndex   %250ms window, 1=active, 0=inactive           %%%
%%%                     fsamp         %sampling freq of EMG                         %%%
%%%                     ISIs          %inter-spike interval (current-prev)          %%%
%%%                     MUPulses      %sample numbers of MU spikes                  %%%
%%%                     IPTs          %estimated spike trains                       %%%
%%%                     PNRs          %pulse-to-noise ratios for each MU            %%%
%%%                     steadyEMGInd     %updated concatenated indeces              %%%
%%%                                      % row1 = 5%, row2 = 20%                    %%%
%%%                     steadyCST     %steady 10s of normalized+smoothed CST        %%%
%%%                                      % row1 = 5%, row2 = 20%                    %%%
%%%                     steadyBinary  %steady 10s of binary MU spike trains         %%%
%%%                                      % row1 = 5%, row2 = 20%                    %%%
%%%                     steadyActiveIndex  %steady 10s of activeIndex               %%%
%%%                                      % row1 = 5%, row2 = 20%                    %%%
%%%                     steadyIPTs    %steady 10s of IPTs                           %%%
%%%                                      % row1 = 5%, row2 = 20%                    %%%
%%%                     steadyPNRs    %PNRs calculated over steady 10s              %%%
%%%                                      % row1 = 5%, row2 = 20%                    %%%
%%%                     SDCST         %SD of CST from 1s epochs of steady 10s       %%%
%%%                                      % row1 = 5%, row2 = 20%                    %%%
%%%                     CVISI         %CV for ISI from 1s epochs of steady 10s      %%%
%%%                                      % row1 = 5%, row2 = 20%                    %%%
%%%                                                                                 %%%
%%%             concatenated          %order = vision, sound, vibration             %%%
%%%                 t5                %target force 5% MVC                          %%%
%%%                     emg           %emg signal trimmed around first ramp         %%%
%%%                     cst           %unedited cumulative spike train              %%%
%%%                     normcst       %normalized to 250ms activeIndex              %%%
%%%                     binary        %binary MU spike trains                       %%%
%%%                     activeIndex   %250ms window, 1=active, 0=inactive           %%%
%%%                     fsamp         %sampling freq of EMG                         %%%
%%%                     ISIs          %inter-spike interval (current-prev)          %%%
%%%                     MUPulses      %sample numbers of MU spikes                  %%%
%%%                     IPTs          %estimated spike trains                       %%%
%%%                     PNRs          %pulse-to-noise ratios for each MU            %%%
%%%                     steadyEMGInd     %updated concatenated indeces              %%%
%%%                                      % row1 = vision                            %%%
%%%                                      % row2 = sound                             %%%
%%%                                      % row3 = vibration                         %%%
%%%                     steadyCST     %steady 10s of normalized+smoothed CST        %%%
%%%                                      % row1 = vision                            %%%
%%%                                      % row2 = sound                             %%%
%%%                                      % row3 = vibration                         %%%
%%%                     steadyBinary  %steady 10s of binary MU spike trains         %%%
%%%                                      % row1 = vision                            %%%
%%%                                      % row2 = sound                             %%%
%%%                                      % row3 = vibration                         %%%
%%%                     steadyActiveIndex  %steady 10s of activeIndex               %%%
%%%                                      % row1 = vision                            %%%
%%%                                      % row2 = sound                             %%%
%%%                                      % row3 = vibration                         %%%
%%%                     steadyIPTs    %steady 10s of IPTs                           %%%
%%%                                      % row1 = vision                            %%%
%%%                                      % row2 = sound                             %%%
%%%                                      % row3 = vibration                         %%%
%%%                     steadyPNRs    %PNRs calculated over steady 10s              %%%
%%%                                      % row1 = vision                            %%%
%%%                                      % row2 = sound                             %%%
%%%                                      % row3 = vibration                         %%%
%%%                     SDCST         %SD of CST from 1s epochs of steady 10s       %%%
%%%                                      % row1 = vision                            %%%
%%%                                      % row2 = sound                             %%%
%%%                                      % row3 = vibration                         %%%
%%%                     CVISI         %CV for ISI from 1s epochs of steady 10s      %%%
%%%                                      % row1 = vision                            %%%
%%%                                      % row2 = sound                             %%%
%%%                                      % row3 = vibration                         %%%
%%%                 t20                                                             %%%
%%%                     ...                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
temp = readcell('E:/Vibration/Participant_data.xlsx','Sheet','General','Range','A2:A169'); %load subject IDs
subjectIDs = strings(length(temp)/6,1);
for i=0:length(subjectIDs)-1
    subjectIDs(i+1,1) = string(temp{i*6+1,1});
end

VIBdata = struct('subject',...
            [struct('data',...
                struct('id',[],...
                'condition',...
                    [struct('t5',[],'t20',[],'t520',[])],...
                'concatenated',...
                    [struct('t5',[],'t20',[])]))]);

for i=1:length(subjectIDs)
    tempStruct = struct('subject',...
        [struct('data',struct(...
            'id',subjectIDs(i),...
            'Vision',struct('t5',struct(...
                                'force',0.0,'emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyForceInd',0.0,'steadyEMGInd',0.0,'EMdelay',0.0,...
                                'steadyForce',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'CVforce',0.0,'SDCST',0.0,'CVISI',0.0),...
                            't20',struct('force',0.0,'emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyForceInd',0.0,'steadyEMGInd',0.0,'EMdelay',0.0,...
                                'steadyForce',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'CVforce',0.0,'SDCST',0.0,'CVISI',0.0),...
                            't520',struct('emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyEMGInd',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'SDCST',0.0,'CVISI',0.0)),...
            'Sound',struct('t5',struct(...
                                'force',0.0,'emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyForceInd',0.0,'steadyEMGInd',0.0,'EMdelay',0.0,...
                                'steadyForce',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'CVforce',0.0,'SDCST',0.0,'CVISI',0.0),...
                            't20',struct('force',0.0,'emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyForceInd',0.0,'steadyEMGInd',0.0,'EMdelay',0.0,...
                                'steadyForce',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'CVforce',0.0,'SDCST',0.0,'CVISI',0.0),...
                            't520',struct('emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyEMGInd',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'SDCST',0.0,'CVISI',0.0)),...
          'Vibration',struct('t5',struct(...
                                'force',0.0,'emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyForceInd',0.0,'steadyEMGInd',0.0,'EMdelay',0.0,...
                                'steadyForce',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'CVforce',0.0,'SDCST',0.0,'CVISI',0.0),...
                            't20',struct('force',0.0,'emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyForceInd',0.0,'steadyEMGInd',0.0,'EMdelay',0.0,...
                                'steadyForce',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'CVforce',0.0,'SDCST',0.0,'CVISI',0.0),...
                            't520',struct('emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyEMGInd',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'SDCST',0.0,'CVISI',0.0)),...
        'concatenated',struct('t5',struct('emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyEMGInd',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'SDCST',0.0,'CVISI',0.0),...
                              't20',struct('emg',0.0,'cst',0.0,'normcst',0.0,...
                                'binary',0.0,'activeIndex',0.0,'fsamp',0.0,...
                                'ISIs',0.0,'MUPulses',0.0,'IPTs',0.0,'PNRs',0.0,...
                                'steadyEMGInd',0.0,'steadyCST',0.0,'steadyBinary',0.0,...
                                'steadyActiveIndex',0.0,'steadyIPTs',0.0,'steadyPNRs',0.0,...
                                'SDCST',0.0,'CVISI',0.0))))...
         ]);
    VIBdata.subject(i) = tempStruct.subject;
end
          