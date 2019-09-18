%% Dots_Delays_preprocess - Filtering - Ciara Devine

clear all
close all



path = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\BDFs\';
matfolder = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\FilteredData\';
Files2 = {'NoLPF_NoHPF','_NoHPF', '_0_05HPF' '_0_15HPF' '_0_25HPF' '_0_35HPF' '_0_45HPF' '_0_55HPF' '_0_65HPF'}; % final part of file name changes depending on which HPF is being used.
SubFold = {'NoLPF_NoHPF\','LPF_NoHPF\','LPF_HPF_0_05\' 'LPF_HPF_0_15\' 'LPF_HPF_0_25\' 'LPF_HPF_0_35\' 'LPF_HPF_0_45\' 'LPF_HPF_0_55\' 'LPF_HPF_0_65\'}; % Which subfolder to look in depending on which high pass filter being used


allsubj=  { 'GD'  'KM'  'OM'  'AF' 'SB' 'DH' 'SM' 'SD' 'HM' 'KW' 'DM' 'AE' 'ST' 'KL' 'SY' 'WR' 'JE' 'CM' 'EC' 'DPM' 'TS' 'SCM'}; % In order of when they were tested
blocks = {[1:5] [1:6] [1:7] [1:10] [1:10] [1:10] [1:10] [1:10] [1:11] [1:9] [1:9] [1:8] [1:8] [1:7] [1:8] [1:9] [1:8] [1:8] [1:8] [1:8] [1:8] [1:8]}; % For HM Block 3 got split into block 3 & 4 due to recording issue

fs=512; % sample rate

chanlocs = readlocs('cap128.loc'); % loads in the channel coordinates, necessary for plotting topographies and for channel interpolation
chanlocs = chanlocs(1:128);

%% Cut off settings for filtering
HPFcutoff=[0 0.05 0.15 0.25 0.35 0.45 0.55 0.65]; % High pass filter cutoff

%% Low Pass
LPFcutoff=35;       % Low Pass Filter cutoff
[Bl, Al] = butter(4,(LPFcutoff*2)/512, 'low');


%% High Pass

HPFcutoff=[0 0.05 0.15 0.25]; % High pass filter cutoff

for Filt=1:length(HPFcutoff)
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
for s=20:length(allsubj); % loop through all subjects
    %% Loop through all blocks
    f = 0; % set file counter to 0
    for b=blocks{s} % loop through all blocks for each session and each subject      
        f = f+1; % update f with each new block
        %% Load in Data
        try
            filename=[path allsubj{s} 'DotsDelays' num2str(b) '-Deci.bdf'];
            EEG = pop_biosig(filename, 'ref',1);
        catch
        end
        try
            filename=[path allsubj{s} 'DotsDelays' num2str(b) '.bdf'];
            EEG = pop_biosig(filename, 'ref',1);
        catch
        end
        try
            filename=[path allsubj{s} 'Dots' num2str(b) '.bdf'];
            EEG = pop_biosig(filename, 'ref',1);
        catch
        end
        try
            filename=[path allsubj{s} 'Dots' num2str(b) '-Deci.bdf'];
            EEG = pop_biosig(filename, 'ref',1);
        catch
        end
        
        %% Resample if sampling rate isn't 512hz 
        if EEG.srate>512 || EEG.srate<512; % srate = sampling rate % Cannot upsample EEG data but this line will produce error if any file is <512., % srate = sampling rate
            EEG = pop_resample(EEG, 512); % if sampling rate is greater than 512 then resample the EEG data at 512.
        end
        
                %%  Detrend in order to remove linear trends in the data
        if detrend_on==1;
            EEG.data=detrend(EEG.data')'; % Detrend works in columns so EEG.data needs to be transposed in order to detrend for each channel. Without transpoing detrend will be done for each time point (across channels) because columns = timepoints as opposed for each channel.
        end
        %save([matfolder SubFold{1} allsubj{s} num2str(b) Files2{1} '_NoLPF'], 'EEG');

        %% Low Pass Filter
        if LPF_on
            EEG.data = filtfilt(Bl,Al,double(EEG.data)')';
        end

        %% Save the data in its current state before looping through filters
        EEGtemp.data = EEG.data;
        
        %% HIGH PASS FILTER
        for Filt = 2:3
            if HPF_on && HPFcutoff(Filt)>0 % filt = 1 or 2 represents no HPF
                EEG.data = filtfilt(Bh{Filt},Ah{Filt},EEGtemp.data')'; 
            end
            
                        %% Save out non CSD data which until this point has been low pass filtered, detrended and potentially high pass filtered
            save([matfolder SubFold{Filt+1} allsubj{s} num2str(b) Files2{Filt+1}], 'EEG');


        end
    end
    
end
