clear;
clc;

folder = uigetdir('Select Folder');
cd(folder);

files = dir('*.mat');
files = struct2cell(files);
files = files(1,:);

fs = 2000;

%%% signal inspection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f = 1:length(files)
    filename = char(files(1,f));
    load(filename);
    L = length(emg_save);
    x1 = emg_save(:,5); x2 = emg_save(:,15); ...
        x3 = emg_save(:,end-10); x4 = emg_save(:,end-5);
    PDS1 = fft(x1); PDS2 = fft(x2);...
        PDS3 = fft(x3); PDS4 = fft(x4);
    P21 = abs(PDS1/L); P22 = abs(PDS2/L);...
        P23 = abs(PDS3/L); P24 = abs(PDS4/L);
    P11 = P21(1:L/2+1); P12 = P22(1:L/2+1);...
        P13 = P23(1:L/2+1); P14 = P24(1:L/2+1);
    P11(2:end-1) = 2*P11(2:end-1); P12(2:end-1) = 2*P12(2:end-1);...
        P13(2:end-1) = 2*P13(2:end-1); P14(2:end-1) = 2*P14(2:end-1);
    freq = fs*(0:(L/2))/L;
    subplot(2,2,1); plot(freq,P11); hold on;...
    subplot(2,2,2); plot(freq,P12); hold on;...
    subplot(2,2,3); plot(freq,P13); hold on;...
    subplot(2,2,4); plot(freq,P14); hold on;
end
%%
%%% signal filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
notch60 = designfilt('bandstopiir','FilterOrder',2,...
            'HalfPowerFrequency1',59,'HalfPowerFrequency2',61,...
            'DesignMethod','butter','SampleRate',fs);
% notch90 = designfilt('bandstopiir','FilterOrder',2,...
%             'HalfPowerFrequency1',89,'HalfPowerFrequency2',91,...
%             'DesignMethod','butter','SampleRate',fs);

for f = 1:length(files)
    filename = char(files(1,f));
    load(filename);
    emg_save = double(emg_save);
    channels = 64;
    emg_notch = zeros(length(emg_save),channels);
    for c = 1:channels
        emg_notch(:,c) = filtfilt(notch60,emg_save(:,c));
        %emg_notch(:,c) = filtfilt(notch90,emg_notch(:,c));
    end
    emg_save = emg_notch;
    %%% smooth edges with window filter %%%
    H = hann(400);  %0.2 second hann window
    leftHann = H(1:200);   %0.1 second left half of hann window
    rightHann = H(201:end);   %0.1 second right half of hann window
    rect = ones(length(emg_save)-length(H),1);   %rect window over EMG duration
    rectHann = [leftHann; rect; rightHann];   %combine hann and rect windows
    for i=1:channels
        emg_save(:,i) = rectHann.*emg_save(:,i);
    end
    
    %%% save %%%
    if exist('Trimmed & Filtered','dir')
    else
        mkdir('Trimmed & Filtered');
    end
    cd('Trimmed & Filtered');
    [~,file,ext] = fileparts(filename);
    filename = sprintf('%s_%s',file,'Filtered');
    save(filename,'emg_save');
    cd(folder);
    clearvars -except f files notch60 folder fs
end
