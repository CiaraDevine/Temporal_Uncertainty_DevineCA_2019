% Difficulties - Behavioural
% Save trigger and latency info from eeg.event into a matfile containing
% simple matrices
% This will make it quick to run the scripts because matlab won't have to
% operate in structural arrays (as in the BDF files). 

clear all 
close all
clc

path = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\BDFs\';
matfolder = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\Events\';


allsubj=  { 'GD'  'KM'  'OM'  'AF' 'SB' 'DH' 'SM' 'SD' 'HM' 'KW' 'DM' 'AE' 'ST' 'KL' 'SY' 'WR' 'JE' 'CM' 'EC' 'DPM' 'TS' 'SCM'}; % In order of when they were tested
blocks = {[1:5] [1:6] [1:7] [1:10] [1:10] [1:10] [1:10] [1:10] [1:11] [1:9] [1:9] [1:8] [1:8] [1:7] [1:8] [1:9] [1:8] [1:8] [1:8] [1:8] [1:8] [1:8]}; % For HM Block 3 got split into block 3 & 4 due to recording issue


for s = 20:length(allsubj); % Scroll through each participant separately 
    
    for b = 1:length(blocks{s}) % For each participant scroll through each block separately (i.e. each file)
        
        clear matfiles trigs trigtimes;
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
        
        if EEG.srate>512 || EEG.srate<512; % srate = sampling rate % Cannot upsample EEG data but this line will produce error if any file is <512.
            EEG = pop_resample(EEG, 512); % if sampling rate is greater than 512 then resample the EEG data at 512.
        end
        matfile{s}{b} = [matfolder, allsubj{s}, 'DotsDelays' num2str(b), '_events.mat']; 
        
        % Extract trigger and latency information from structure format to
        % matrix format 
        
        for i=1:length(EEG.event)
        trigs(i)=EEG.event(i).type;  % pull out trigger labels and, below, trigger times from the EEGlab file structure. 
        % This step could be skipped but is slightly handier than having to call on the specific subdirectories each time
        trigtimes(i)=round(EEG.event(i).latency);
        end
        
 save(matfile{s}{b}, 'trigs' , 'trigtimes')  % Each time matlab gets to the end of 
 % a block loop it will save a matfile for that block for the given
 % subject.The saved variables and files must then be cleared at the
 % beginning of the next block loop. 
 
    end

end 
        
 
