clear;
clc;

trimFolder = 'E:\Vibration\To Trim';
cd(trimFolder);
subFolders = dir;
subFolders = subFolders(3:length(subFolders));

for s = 1:length(subFolders)
    cd(subFolders(s).name);
    files = dir('*.mat');
    files = struct2cell(files);
    files = files(1,:);
    for f=1:length(files)
        filename = char(files(1,f));
        load(filename);
        
        plot(SIG{5,3});
        title(filename);
        [trim,y] = ginput(1);
        trim = round(trim);
        SIGlength = trim/2000;
        stopSIGInt = SIGlength;
        IPTs = IPTs(:,1:trim);
        
        H = hann(400);  %0.2 second hann window
        leftHann = H(1:200);   %0.1 second left half of hann window
        rightHann = H(201:end);   %0.1 second right half of hann window
        rect = ones(trim-length(H),1);   %rect window over EMG duration
        rectHann = [leftHann; rect; rightHann];   %combine hann and rect windows
        
        for r=1:13
            for c=1:5
                if ~isempty(SIG{r,c})
                    SIG{r,c} = SIG{r,c}(1,1:trim);
                    SIG{r,c} = rectHann'.*SIG{r,c};
                end
            end
        end

        for mu = 1:numel(MUPulses)
            for i = 1:size(MUPulses{1,mu},2)
                if MUPulses{1,mu}(1,i)>trim
                    if i==1
                        MUPulses{1,mu}=[];
                        break;
                    else
                        MUPulses{1,mu} = MUPulses{1,mu}(1:i-1);
                        break;
                    end
                end
            end
        end
        [~,file,ext] = fileparts(filename);
        filename = sprintf('%s_%s',file,'FirstRamp');
        save(filename);
        clearvars -except s f files trimFolder subFolders
    end
    cd(trimFolder);
end