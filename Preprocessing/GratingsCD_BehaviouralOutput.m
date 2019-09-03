% Behavioural Script for training study - Ciara Devine

% Task - two alternative forced choice contrast discrimination task - superimposed grating stimuli.
% Subjects required to detect change in the contrast of the gratings (one grating increases in contrast
% while the other decreases). Subjects are required to indicate the direction of the emerging target
% as quickly as possible once detected (ie leftward or rightward) by
% clicking either the left or right mouse button with the corresponding
% thumb.

clear all
close all
clc

Datafolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\Events\'];

subjID =  { 'CJG'  'CR'  'AR'  'KD'  'TB'  'AOB'  'PM'  'CMG'  'SR'  'AB' 'GK' 'ED' 'SC' 'CB' 'CE' 'SH' 'RMC' 'JH', 'JR', 'ROC' 'KM' 'ROS'}; % In order of when they were tested 
sessions = {[1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:4] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5], [1:5] [1:5] [1:5] [1:5]}; 
blocks =  { {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:29; 30:39; 40:49}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50},  {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50}};
allblocks_sub={[50] [49] [50] [50] [50] [50] [50] [50] [50] [50] [40] [50] [50] [50] [50] [50] [50] [50] [50] [50] [50] [50]};
subjects=22;
allblocks = 50;
blocks_persess = 10; % 10 blocks in seach session
allsessions = 5;
Max_trials_persub = 500*5; % Each subject did completely a maximum of 2500 trials. 

targs = [8 9]; % 8 = left targ, 9 = right targ
resps = [12 13]; % 12 = left resp, 13 = right resp
rewards(1,150:2001)=80:-0.0216:40; %% Rewards matrix goes from 1 - 2000ms. 0 points given for RT up to 150ms (ie RT<150 not accepted as correct).
% From RT = 150ms onwards points drop from 80 to 40 across target duration in steps of .0216. Why .0216??? ...
% 2000ms-150ms = 1850ms. 80-40 = 40. 1850/40 = .0216. Thus points are reduced by .0216 with every ms because reducing by this value brings reward from 80 to 40 points over the course of the 1850ms.


%% 
nn= 0; % trial counter for full dataset 
% create this outside of the subjects loop so that it is not over written with each new subject (or session or block). 
% As each new trial in processed this will be updated to reflect the number of trials that have come before 
for s = 1:length(subjID); % Loop through each subject
    
    SubTrialNo = 0; % C
    
    for d = 1:length(sessions{s}) % Loop through each session for each subject
        
        
        f=0; % set files counter to 0 for each new session. Files counter serves as a within session block counter
        
        for b = blocks{s}{d} % Loop through each block (i.e. each events file) for each session and subject
            
            % Using {} tells matlab to go inside s and inside d (i.e. to pull out all that is inside d (ie the blocks))

            f = f+1; % Within each session f will update in accordance with each new block (eg block 11 will have f=1 because its the first file/block of the second session)
            
            filename = [Datafolder subjID{s} num2str(b) '_events.mat'];
            files{s}{d}{f} = [Datafolder subjID{s} num2str(b) '_events.mat']; % create a cell arrary containing all file names - sanity check to ensure block loop is working as it should.
            load(filename); % Load in the file corresponding to the current block ie block{s}{d}
            


            for i = 1:length(trigs) % loop through each trigger in the file. Trigs is saved in the events file just loaded
                early=0; 
                try 
                    
                if ismember(trigs(i), [8 9]) && ismember(trigs(i+1),[12 13]); % Check for & process trials where the response was made after the target 
                    SubTrialNo = SubTrialNo+1; % Count number of trials this subject now has  
                    nn = nn +1; % Update the number of trials 
                    SubTrials(nn) = SubTrialNo; 
                    SubID(nn)=s; % Record what subject this trial belongs to
                    RT(nn) = round((trigtimes(i + 1) - trigtimes(i))*1000/512); % To convert from samples to ms you must multiply by 1000/512 % 512 is the sampling rate
                    Sessid(nn) = d;  % record what session this trial belongs to               
                    Block(1,nn) = f; % record what block within a session this trial belongs to 
                    Block(2,nn) = b; % record what bloack this trial belongs to ignoring session 
                    Targ(nn) = trigs(i); % Record what target was presented on this block
                    Resp(nn) = trigs(i+1); % Record what response was made on this trial 
                       try Delay(1,nn) = trigtimes(i)-trigtimes(i-1); % Calculate how long (in samples) the delay was on this trial and then below, categorise the delays according to whether they were ~800/1200/1600ms
                    if Delay(1,nn)<500, Delay(2,nn)=1; elseif Delay(1,nn)>500 && Delay(1,nn)<700, Delay(2,nn)=2;elseif Delay(1,nn)>700 && Delay(1,nn)<900, Delay(2,nn)=3; end; catch; end;
                
                if (ismember(trigs(i), [8]) && ismember(trigs(i+1),[12])) || (ismember(trigs(i), [9]) && ismember(trigs(i+1),[13])); 
                    Accuracy(nn) = 1; % If a left target was followed by a left response or vice versa for right targets, mark this trial as correct
                    Score(nn) = round(rewards(RT(nn))); % Calculate the score for this trial based on RT 
                elseif (ismember(trigs(i), [8]) && ismember(trigs(i+1),[13])) || (ismember(trigs(i), [9]) && ismember(trigs(i+1),[12])); 
                    Accuracy(nn) = 0; % If the incorrect response was given mark this trial as incorrect using a 0 
                    Score(nn) =0; % No points for incorrect trials 
                end
                
                elseif ismember(trigs(i), [12 13]) && ismember(trigs(i+1),[8 9]); % check for trials in which the response came before the target - Process these as false alarms
                    %
                    early=1; 
                    SubTrialNo = SubTrialNo+1; 
                    nn = nn +1; % Update the number of trials 
                    SubTrials(nn) = SubTrialNo; 
                    SubID(nn)=s; 
                    RT(nn) = round((trigtimes(i) - trigtimes(i+1))*1000/512);
                    if RT(nn)<-600, 
                        stop=1; 
                    end
                    Sessid(nn) = d;                 
                    Block(1,nn) = f;
                    Block(2,nn) = b;
                    Targ(nn) = trigs(i+1);
                    Resp(nn) = trigs(i);
                    Accuracy(nn) = 3; 
                     Score(nn) =0;
                       try Delay(1,nn) = trigtimes(i+1)-trigtimes(i-1); 
                    if Delay(1,nn)<500, Delay(2,nn)=1; elseif Delay(1,nn)>500 && Delay(1,nn)<700, Delay(2,nn)=2;elseif Delay(1,nn)>700 && Delay(1,nn)<900, Delay(2,nn)=3; end; catch; end;
                
                
                elseif (trigs(i) == 8 && trigs(i+1) == 5) || (trigs(i) == 9 && trigs(i+1) == 5) && ~early % Check for trials where no response was made 
                    
                    SubTrialNo = SubTrialNo+1; 
                    nn = nn +1; % Update the number of trials 
                    SubTrials(nn) = SubTrialNo; 
                    SubID(nn)=s; 
                    Accuracy(nn) = 2;
                    RT(nn) = NaN;
                    Score(nn) = 0;
                    Sessid(nn) = d;
                    Block(1,nn) = f;
                    Block(2,nn) = b;
                    Resp(nn) = NaN;
                    try Delay(1,nn) = trigtimes(i)-trigtimes(i-1); 
                    if Delay(1,nn)<500, Delay(2,nn)=1; elseif Delay(1,nn)>500 & Delay(1,nn)<700, Delay(2,nn)=2; elseif Delay(1,nn)>700 & Delay(1,nn)<900, Delay(2,nn)=3; end; catch; end;
                    SubID(nn) = s;
                    Targ(nn) = trigs(i);
                    
                end   
                catch 
                end
                
            end

% Hack Job to deal with known abnormalities in this data 
% Remove Score values for blocks where 50 trials were not completed 
Score(SubID==2 & Block(2,:)==9)=NaN;
Score(SubID==2 & Block(2,:)==21) = NaN;
Score(SubID==2 & Block(2,:)==22) = NaN; 
Score(SubID==1 & Block(2,:)==1) = NaN; 

                
        end


    end
    %  Organise data in format suitable for later modelling 
conds =[1:3; 4:6]; 
Data = [Delay(2,SubID==s).*(Targ(SubID==s)); Accuracy(SubID==s); RT(SubID==s)]; 
Data(1,Data(1,:)==8)=1; Data(1,Data(1,:)==16)=2; Data(1,Data(1,:)==24)=3; 
Data(1,Data(1,:)==9)=4; Data(1,Data(1,:)==18)=5; Data(1,Data(1,:)==27)=6; 
% save subject specific data
save([Datafolder num2str(subjID{s}) 'ModellingData'],'Data');

end
% 
RT(RT<-600)=NaN;
% Save output 
save([Datafolder 'BehaviouralOutput'],'SubID', 'Accuracy', 'RT', 'Score', 'Sessid', 'Block', 'Delay','Targ', 'Resp', 'SubTrials');

% Save full dataset 
Data=[Delay(2,:).*(Targ); Accuracy; RT];
Data(1,Data(1,:)==8)=1; Data(1,Data(1,:)==16)=2; Data(1,Data(1,:)==24)=3; 
Data(1,Data(1,:)==9)=4; Data(1,Data(1,:)==18)=5;Data(1,Data(1,:)==27)=6; 
save([Datafolder 'ModellingData'],'Data');  

           
            