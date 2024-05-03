 %%% clear old variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc;

%%% import OTB matlab file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[FileName,PathName] = uigetfile('*.mat','Select mat file');
data_struct = load(fullfile(PathName,FileName));
data = cell2mat(data_struct.Data); % last column is AUX trigger
time = cell2mat(data_struct.Time);
numChannels = size(data,2)-1;

%%% plot AUX trigger channel to find trim_start time %%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1)
plot(time,data(:,numChannels+1))
legend('trigger')
subplot(2,1,2)
plot(time,data(:,15))
legend('EMG')
datacursormode on
%%
%%% trim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trim_start = 5.4175; % update this value with trim start time
trim_end = 32.9345; %time(end); % update this value with trim end time

% find array indeces to be trimmed
i_trim_start = 1; % array index of trim start
for i=1:length(time)
    if abs(time(i) - trim_start) < 1e-04   % boolean error workaround
        i_trim_start = i;
        break
    end
end

i_trim_end = length(time); % array index of trim end
for i=1:length(time)
    if abs(time(i) - trim_end) < 1e-04   % boolean error workaround
        i_trim_end = i;
        break
    end
end

% populate trimmed arrays
emg_trim = zeros(length(time(i_trim_start:i_trim_end)),numChannels);
time_trim = time(i_trim_start:i_trim_end);
for i=1:numChannels
    emg_trim(:,i) = data(i_trim_start:i_trim_end,i);
end

%%% smooth edges with window filter and remove DC offset %%%%%%%%%%%%%%%%%%
emg_windowed = zeros(size(emg_trim));
H = hann(400);  %0.2 second hann window
leftHann = H(1:200);   %0.1 second left half of hann window
rightHann = H(201:end);   %0.1 second right half of hann window
rect = ones(length(time_trim)-length(H),1);   %rect window over EMG duration
rectHann = [leftHann; rect; rightHann];   %combine hann and rect windows

for i=1:numChannels
    emg_trim(:,i) = emg_trim(:,i)-mean(emg_trim(:,i));
    emg_windowed(:,i) = rectHann.*emg_trim(:,i);
end

%%% save trimmed emg signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[file,path] = uiputfile('*.mat','Save file name');
fullFileName = fullfile(path,file);
emg_save = emg_windowed;
save(fullFileName,'emg_save');
