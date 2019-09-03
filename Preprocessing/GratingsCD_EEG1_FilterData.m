%% Gratings Training Delays_preprocess - Filtering - Ciara Devine

clear all
close all



path = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\GratingsCD_Data\RawDataBDFs\'];
matfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\FilteredData\']; % Where to look for incoming data
Files2 = {'_NoHPF', '_0_05HPF' '_0_15HPF' '_0_25HPF' '_0_35HPF' '_0_45HPF' '_0_55HPF' '_0_65HPF'}; % final part of file name changes depending on which HPF is being used.
SubFold = {'LPF_NoHPF\','LPF_HPF_0_05\' 'LPF_HPF_0_15\' 'LPF_HPF_0_25\' 'LPF_HPF_0_35\' 'LPF_HPF_0_45\' 'LPF_HPF_0_55\' 'LPF_HPF_0_65\'}; % Which subfolder to look in depending on which high pass filter being used


allsubj =  { 'CJG'  'CR'  'AR' 'KD'  'TB'  'AOB'  'PM'  'CMG'  'SR'  'AB' 'GK' 'ED' 'SC' 'CB' 'CE' 'SH' 'RMC' 'JH' 'JR' 'ROC' 'KM' 'ROS'}; % In order of when they were tested; KD and JR: too many trials lost due to blinks - excluded!
sessions = {[1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:4] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5]};
blocks =  { {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:29; 30:39; 40:49}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50} ,{1:10; 11:20; 21:30; 31:40; 41:50}};

fs=512; % sample rate

chanlocs = readlocs('cap128.loc'); % loads in the channel coordinates, necessary for plotting topographies and for channel interpolation
chanlocs = chanlocs(1:128);

%%         Make Filters

%% Cut off settings for filtering
LPFcutoff=35;       % Low Pass Filter cutoff
HPFcutoff=[0 0.05 0.15 0.25 0.35 0.45 0.55 0.65]; % High pass filter cutoff

%% Low Pass
[Bl, Al] = butter(2,(LPFcutoff*2)/fs, 'low'); % (LPFcutoff*2)/fs is equivalent mathematically to LPFcutoff/(fs/2) where fs/2=nyquist


%% High Pass
HPFStrengths=[2]; % Indices for the cut offs above. 

for Filt=HPFStrengths
    clear bh ah
    if HPFcutoff(Filt)>0 %Design high pass filter
        [bh,ah] = butter(2,(HPFcutoff(Filt)*2)/fs,'high');
        Bh{Filt}=bh; Ah{Filt}=ah;
    else Bh{Filt}=NaN; Ah{Filt}=NaN;
    end
end

%% Turning on/off filters
LPF_on = 1;    % To turn on LPF set LPF = 1. To turn off LPF set LPF = 0.
detrend_on = 1; % To turn on detrending function set detrend = 1. to turn off, set to 0.
HPF_on = 1;



%% Loop Through All Subjects
for s=1:length(allsubj); % loop through all subjects
        for d=1:5 
            %% Loop through all blocks
            f = 0; % set file counter to 0
            for b=blocks{s}{d} % loop through all blocks for each session and each subject
                f = f+1; % update f with each new block
                filename=[path allsubj{s} num2str(b) '.bdf'];
                EEG = pop_biosig(filename, 'ref',1);
                
                if EEG.srate>512 || EEG.srate<512; % srate = sampling rate % Cannot upsample EEG data but this line will produce error if any file is <512., % srate = sampling rate
                    EEG = pop_resample(EEG, 512); % if sampling rate is greater than 512 then resample the EEG data at 512.
                end
                
                %% Low Pass Filter
                if LPF_on
                    EEG.data = single(filtfilt(Bl,Al,double(EEG.data)')');
                end
                
                %%  Detrend in order to remove linear trends in the data
                if detrend_on==1;
                    EEG.data=detrend(EEG.data')'; % Detrend works in columns so EEG.data needs to be transposed in order to detrend for each channel. Without transpoing detrend will be done for each time point (across channels) because columns = timepoints as opposed for each channel.
                end
                
                %% Save the data in its current state before looping through filters
                EEGtemp.data = EEG.data;
                
                %% HIGH PASS FILTER
                for Filt = HPFStrengths
                    if HPF_on && HPFcutoff(Filt)>0 % filt = 1 or 2 represents no HPF
                        EEG.data = single(filtfilt(Bh{Filt},Ah{Filt},double(EEGtemp.data)')');
                    end
                    
                    
                    %% Save out non CSD data which until this point has been low pass filtered, detrended and potentially high pass filtered
                    save([matfolder '\' SubFold{Filt} allsubj{s} num2str(b) Files2{Filt}], 'EEG');
                    % Filt = 1 = Entirely Unfiltered data - saved out above
                    % already before LPF and detrend
                    
                end
            end
        end
end
