clear all
close all
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% OTBio EMG Signal %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 2000;
NumChannels = 64;
%%% Nov 21 grid setup %%%
% electrode_map = [25 26 27 28 29 30 31 32; % distal grid
%                  17 18 19 20 21 22 23 24;
%                   9 10 11 12 13 14 15 16;
%                   1  2  3  4  5  6  7  8;
%                  33 34 35 36 37 38 39 40; % proximal grid
%                  41 42 43 44 45 46 47 48;
%                  49 50 51 52 53 54 55 56;
%                  57 58 59 60 61 62 63 64];
%%% Dec 19 grid setup %%%
electrode_map = [25 26 27 28 29 30 31 32; % distal grid
                 17 18 19 20 21 22 23 24;
                  9 10 11 12 13 14 15 16;
                  1  2  3  4  5  6  7  8;
                 57 58 59 60 61 62 63 64; % proximal grid
                 49 50 51 52 53 54 55 56;
                 41 42 43 44 45 46 47 48;
                 33 34 35 36 37 38 39 40];
x_coord = [1 2 3 4 5 6 7 8];
y_coord = [1 2 3 4 5 6 7 8];

%%% import OTB matlab file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[file_name,path_name] = uigetfile('*.mat','Select mat file');
data_struct = load(fullfile(path_name,file_name));
data = double(cell2mat(data_struct.Data));
time = double(cell2mat(data_struct.Time));

% eliminate DC offset from signals
for i=1:NumChannels
data(:,i) = data(:,i)-mean(data(:,i));
end

% bandpass filter signals
[b,a] = butter(4,[20 500]/(fs/2));
data = (filtfilt(b,a,data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Amplitude Maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,file,ext] = fileparts(file_name);
file_name = sprintf('%s_%s',file,'AmplitudeMaps');
v = VideoWriter(file_name);
v.Quality = 100;
v.FrameRate = 8;
open(v);
epoch_length = 0.5;  % length of the epoch for RMS calculation (seconds)
final_epoch = floor(size(data,1)/fs);
max_RMS = 0;
fg = figure(1);

for e = 0:epoch_length:final_epoch-1
    epoch_indeces = [(e*fs)+1:(e+epoch_length)*fs+1];
    % computes the RMS over the selected epoch
    for row = 1:size(electrode_map,1)
        for col = 1:size(electrode_map,2)
            RMS_map(row,col) = sqrt(mean(data(epoch_indeces,electrode_map(row,col)).^2));
            if RMS_map(row,col) > max_RMS
                max_RMS = RMS_map(row,col);
            end
        end
    end
    [bcX, bcY] = barycenter2D(x_coord,y_coord,RMS_map);
    
    % plot amplitude map
    set(0,'CurrentFigure',fg);
    imagesc(RMS_map,[0 800])
    colormap(jet)
    colorbar
    hold on
    plot(bcX,bcY,'+w','MarkerSize',20,'LineWidth',2)
    hold off
    title(['RMS map from T = ' num2str(e) 's to T = ' num2str(e+epoch_length) 's'])
    xlabel('Columns')
    ylabel('Rows')
    set(gca,'xtick',[1:size(electrode_map,2)])
    set(gca,'ytick',[1:size(electrode_map,1)])
    frame = getframe(gcf);
    writeVideo(v,frame);

%    % plot interpolated amplitude map 
%     interp_factor = 2; % interpolation factor
%     RMS_map_int = interp2(RMS_map,interp_factor);
%     set(0,'CurrentFigure',fg);
%     colormap(jet)
%     colorbar
%     title(['Interpolated RMS map from T = ' num2str(e) 's to T = ' num2str(e+epoch_length) 's'])
%     xlabel('Columns')
%     ylabel('Rows')
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end
close(fg);
close(v);

function [bcX, bcY] = barycenter2D(x,y,C)
%% BARYCENTER2D returns x and y barycenter values of density map C
%  Assumption: Input coordinated x and y are evenly spaced!
%
%   INPUT (from Matlab(R))
%     x    ... vector of length n, x-coordinates of image
%     y    ... vector of length m, y-coordinates of image
%     C    ... matrix of size n x m, density value for each coordiante
%
%   OUTPUT
%     bcX  ... double, x-coordinate value of the barycenter
%     bcY  ... double, y-coordinate value of the barycenter
%
  % Find the center of mass of matrix C // [2]
    C = C/sum(C(:));
    [m,n]=size(C);
    [I,J]=ndgrid(1:m,1:n);
    iX = dot(J(:),C(:));
    iY = dot(I(:),C(:));
  % Transform the index of the weights to coordinates
    xIdx = 1:size(C,2); % x axis index vector
    yIdx = 1:size(C,1); % y axis index vector
  % align true x axis with index vector
    pX(1) = sum((xIdx-mean(xIdx)).*(x-mean(x)))./sum((xIdx-mean(xIdx)).^2); % slope
    pX(2) = mean(x)-pX(1)*mean(xIdx); % intercept
    bcX = pX(1)*iX+pX(2); % eval barycenter x with slope and intercept
  % align true y axis with index vector
    pY(1) = sum((yIdx-mean(yIdx)).*(y-mean(y)))./sum((yIdx-mean(yIdx)).^2); % slope
    pY(2) = mean(y)-pY(1)*mean(yIdx); % intercept
    bcY = pY(1)*iY+pY(2); % eval barycenter x with slope and intercept
end
