clear all 
close all
clc
eeglab; 
Path = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\GratingsCD_Data\RawDataBDFs\'];
Eventsfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\Events\']; 

subjID = {}; %{'CJG' 'CR' 'AR' 'KD' 'TB' 'AOB' 'PM' 'SR' 'AB' 'GK' 'ED' 'SC' 'CB' 'CE''SH' 'RMC''JH''JR' 'ROC' 'KM' 'ROS'};
allblocks = { }; % {[1:50] [1:49] [1:50] [1:50] [1:50] [1:50] [1:50] [1:50] [1:50] [1:50] [1:50] [1:50] [1:50] 1:50] [1:50] [1:50] [1:50]  [1:50]  [1:50]};

for s = 1:length(subjID); % Scroll through each participant separately 
    
    for b = 1:length(allblocks{s}) % For each participant scroll through each block separately (i.e. each file)
        
        clear matfiles trigs trigtimes;
        filename = [subjID{s} num2str(b) '.bdf'];
        
        EEG = pop_biosig(filename); 
        
        matfile{s}{b} = [Eventsfolder, subjID{s}, num2str(b), '_events.mat']; 
        
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
        
 
