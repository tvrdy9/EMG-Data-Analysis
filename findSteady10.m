clearvars -except VIBdata;
clc;

fc = 20;         % Filter cutoff = 20 Hz
fsamp = 2000;           % Sampling frequency (rate) = 2000 Hz
[filt2, filt1] = butter(4,fc/(fsamp/2));      % Apply Butterworth filter
TrialNum = double(readmatrix('E:/Vibration/Participant_data.xlsx','Sheet','General','Range','D2:D169'));
f=87;

for s = 15 %1:length(VIBdata.subject)
    for c = 2 %1:3
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
        for t=5 %[5,20]
            %%% Find Steady 10s %%%
            % load Force file
            targetForce = append('t',string(t));
            rawforce = VIBdata.subject(s).data.(conditionString).(targetForce).force;
            rawtime = 0:1/2000:length(rawforce)-1/2000;
            rawtime = rawtime(end);

            % filter
            force = filtfilt(filt2, filt1, rawforce);

            % plot
            forcefig1 = figure(1);
            forceplot1 = plot(force,'LineWidth',2,'Color',[0.6 0.6 0.6]);
            set(forcefig1, 'Position', [400, 250, 1200, 400])
            ylabel('\fontsize{12}Force (Volts)');
            xlabel('\fontsize{12}Time');
            title('\fontsize{12}Select values of interest: Baseline, Target:')

            % User Input - Select start & end of trial & check that the section is long enough
            [baseline1a,~] = ginput(1);
            [baseline1b,~] = ginput (1);  

            [start1,~] = ginput(1);
            [end1,~] = ginput(1);

            hold on  

            % Adjust start & end points to be integers
            start1 = floor(start1);
            end1 = floor(end1);

            % Steadiest window parameters for 10-s window
            time = 10;      % Number of seconds desired for window (10s)
            window = time*fsamp; % Convert time to # data points
            jump = fsamp/time;   % How far apart each 10s window should be (0.1s)

            % Preallocate vectors
            mean1 = zeros(10000,1);
            sd1 = zeros(10000,1);
            cv1 = zeros(10000,1);

            % Variables in for loop
            r1 = 1; % Counter for Target phase
            section1 = force(start1+fsamp:start1+fsamp+window-1); % Vector with first 10s window
            sweeps1 = floor(((end1-2*fsamp-start1)-window)/jump); % Number of loops needed

            % Define baseline
            baseline1=force(floor(baseline1a):floor(baseline1b));

            % For Loop to scan for steadiest 10s window
            for j = 1:sweeps1
                tempmean1 = mean(section1)-mean(baseline1);
                mean1(r1,1) = tempmean1;
                tempsd1 = std(section1);
                sd1(r1,1) = tempsd1;
                cv1(r1,1) = coefvar(tempmean1,tempsd1);
                r1 = r1 + 1;
                section1 = force(start1+fsamp+(jump*j):start1+fsamp+(jump*j)+window-1);
            end

            % Find lowest SD & CV
            % Trim zeros off end of vectors
            trimsd1 = find(sd1,1,'last');
            sd1 = sd1(1:trimsd1);
            trimcv1 = find(cv1,1,'last');
            cv1 = cv1(1:trimcv1);

            % Minimums  (SN = sweep number)   
            [minsd1,SNsd1] = min(sd1);
            [mincv1,SNcv1] = min(cv1);

            % Plot results
            % Selected portion of force trace with minimum sd and cv (Should pretty much always be the same for sd and CV!)
            sdIndex1 = (start1+fsamp+SNsd1*jump);
            selectedsd1 = force(sdIndex1:sdIndex1+window-1);
            cvIndex1 = (start1+fsamp+SNcv1*jump);
            selectedcv1 = force(cvIndex1:cvIndex1+window-1);

            % Align selected portion to x-axis
            align1sd = (sdIndex1:sdIndex1+window-1)';
            plot1sd = horzcat(align1sd,selectedsd1);

            % Align selected portion to x-axis
            align1cv = (cvIndex1:cvIndex1+window-1)';
            plot1cv = horzcat(align1cv,selectedcv1);

            % Mean force during this window
            minmean1 = mean1(SNsd1);

            % Plot
            forcefig1 = figure(1);
            forceplot1 = plot(force,'LineWidth',2,'Color',[0.6 0.6 0.6]);  % Original
            hold on

            plot(plot1cv(:,1),plot1cv(:,2),'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
            y1 = get(gca,'ylim');
            plot([cvIndex1, cvIndex1],y1,'LineWidth',0.8,'Color',[0.8500 0.3250 0.0980],'LineStyle','--')
            plot([cvIndex1+window-1, cvIndex1+window-1],y1,'LineWidth',0.8,'Color',[0.8500 0.3250 0.0980],'LineStyle','--')

            title('\fontsize{14}Force Steadiness')
            hold off  

            % Start and endtime
            startSamp = cvIndex1;
            startTime = startSamp/fsamp;
            
            %%% Save Steady Indeces to Excel %%%
            sampCellString = sprintf('A%d',f+1);
            timeCellString = sprintf('B%d',f+1);
            forceCellString = sprintf('C%d',f+1);
            sdCellString = sprintf('D%d',f+1);
            cvCellString = sprintf('E%d',f+1);
            filename = "test.xlsx";
            pathname = "E:/Vibration";
            filename = fullfile(pathname,filename);
            writematrix(startSamp,filename,'Sheet','General','Range',sampCellString);
            writematrix(startTime,filename,'Sheet','General','Range',timeCellString);
            writematrix(minmean1,filename,'Sheet','General','Range',forceCellString);
            writematrix(minsd1,filename,'Sheet','General','Range',sdCellString);
            writematrix(mincv1,filename,'Sheet','General','Range',cvCellString);
            f=f+1;
        end
    end
end