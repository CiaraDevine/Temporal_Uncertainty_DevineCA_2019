%% Behavioural Script for Dots Study - Ciara Devine


clear all
close all

Datafolder = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\Events\';


subjID=  { 'GD'  'KM'  'OM'  'AF' 'SB' 'DH' 'SM' 'SD' 'HM' 'KW' 'DM' 'AE' 'ST' 'KL' 'SY' 'WR' 'JE' 'CM' 'EC' 'DPM' 'TS' 'SCM'}; % In order of when they were tested
blocks = {[1:5] [1:6] [1:7] [1:10] [1:10] [1:10] [1:10] [1:10] [1:11] [1:9] [1:9] [1:8] [1:8] [1:7] [1:8] [1:9] [1:8] [1:8] [1:8] [1:8] [1:8] [1:8]}; % For HM Block 3 got split into block 3 & 4 due to recording issue


rewards(1,150:2001)=80:-0.0216:40; %% Rewards matrix goes from 1 - 2000ms. 0 points given for RT up to 150ms (ie RT<150 not accepted as correct).
% From RT = 150ms onwards points drop from 80 to 40 across target duration in steps of .0216. Why .0216??? ...
% 2000ms-150ms = 1850ms. 80-40 = 40. 1850/40 = .0216. Thus points are reduced by .0216 with every ms because reducing by this value brings reward from 80 to 40 points over the course of the 1850ms.

allsubs=[1:22];
duds=[1:3];
subs= ~ismember(allsubs, duds);

%% 
AllTargs= 0; % create this outside of the subjects loop so that it is not over written with each new subject (or session or block). The idea is to keep track of ALL trials, from all subjects
for s = 1:length(subjID); % Loop through each subject
    
    SubTrialNo = 0; Data1=[]; 
    
        
        for b = 1:length(blocks{s}) % Loop through each block (i.e. each events file) for each session and subject
            
            % Using {} tells matlab to go inside s and inside d (i.e. to pull out all that is inside d (ie the blocks))

            
            filename = [Datafolder subjID{s} 'DotsDelays'  num2str(b) '_events.mat'];
            files{s}{b} = [Datafolder subjID{s} 'DotsDelays'  num2str(b) '_events.mat']; % create a cell arrary containing all file names - sanity check to ensure block loop is working as it should.
            load(filename); % Load in the file corresponding to the current block ie block{s}{d}
            


            for i = 1:length(trigs) % loop through each trigger in the file. Trigs is saved in the events file just loaded
                early=0; 
                if ismember(trigs(i), [8 9 12 13 16 17]) && ismember(trigs(i+1),[20 25]) % && round((trigtimes(i + 1) - trigtimes(i))*1000/512) <= 2000 %Only include trials where RT was within 2 seconds - eliminate any trials where a blip in task may have caused the length of the target to be longer than 2seconds.
                    SubTrialNo = SubTrialNo+1; 
                    AllTargs = AllTargs +1; % Update the number of trials 
                    SubTrials(AllTargs) = SubTrialNo; 
                    SubID(AllTargs)=s; 
                    RT(AllTargs) = round((trigtimes(i + 1) - trigtimes(i))*1000/512);
                    Block(2,AllTargs) = b;
                    Targ(AllTargs) = trigs(i);
                    Resp(AllTargs) = trigs(i+1);
                    try Delay(1,AllTargs) = trigtimes(i)-trigtimes(i-1); 
                    if Delay(1,AllTargs)<500, Delay(2,AllTargs)=1; elseif Delay(1,AllTargs)>500 && Delay(1,AllTargs)<700, Delay(2,AllTargs)=2;elseif Delay(1,AllTargs)>700 && Delay(1,AllTargs)<900, Delay(2,AllTargs)=3; end; catch; end;
                
                if (ismember(trigs(i), [8 12 16]) && ismember(trigs(i+1),[20])) || (ismember(trigs(i), [9 13 17]) && ismember(trigs(i+1),[25])); 
                    Accuracy(AllTargs) = 1; % Mark this trial as correct
                    Score(AllTargs) = round(rewards(RT(AllTargs)));
                elseif (ismember(trigs(i), [8 12 16]) && ismember(trigs(i+1),[25])) || (ismember(trigs(i), [9 13 17]) && ismember(trigs(i+1),[20])); 
                    Accuracy(AllTargs) = 0; % Mark this trial as correct
                    Score(AllTargs) =0;
                end
                
                elseif ismember(trigs(i), [20 25]) && ismember(trigs(i+1),[ 8 9 12 13 16 17])% && round((trigtimes(i + 1) - trigtimes(i))*1000/512) <= 2000 %Only include trials where RT was within 2 seconds - eliminate any trials where a blip in task may have caused the length of the target to be longer than 2seconds.
                    early=1; 
                    SubTrialNo = SubTrialNo+1; 
                    AllTargs = AllTargs +1; % Update the number of trials 
                    SubTrials(AllTargs) = SubTrialNo; 
                    SubID(AllTargs)=s; 
                    RT(AllTargs) = round((trigtimes(i) - trigtimes(i+1))*1000/512);
                    Block(2,AllTargs) = b;
                    Targ(AllTargs) = trigs(i+1);
                    Resp(AllTargs) = trigs(i);
                    Accuracy(AllTargs) = 3; 
                    Score(AllTargs) =0;
                       try Delay(1,AllTargs) = trigtimes(i+1)-trigtimes(i-1); 
                    if Delay(1,AllTargs)<500, Delay(2,AllTargs)=1; elseif Delay(1,AllTargs)>500 && Delay(1,AllTargs)<700, Delay(2,AllTargs)=2;elseif Delay(1,AllTargs)>700 && Delay(1,AllTargs)<900, Delay(2,AllTargs)=3; end; catch; end;
                
                
                elseif ismember(trigs(i), [8 9 12 13 16 17]) && trigs(i+1) == 5 && ~early
                    
                    SubTrialNo = SubTrialNo+1; 
                    AllTargs = AllTargs +1; % Update the number of trials 
                    SubTrials(AllTargs) = SubTrialNo; 
                    SubID(AllTargs)=s; 
                    Accuracy(AllTargs) = 2;
                    RT(AllTargs) = NaN;
                    Score(AllTargs) = 0;
                    Block(2,AllTargs) = b;
                    Resp(AllTargs) = NaN;
                    try Delay(1,AllTargs) = trigtimes(i)-trigtimes(i-1); 
                    if Delay(1,AllTargs)<500, Delay(2,AllTargs)=1; elseif Delay(1,AllTargs)>500 & Delay(1,AllTargs)<700, Delay(2,AllTargs)=2; elseif Delay(1,AllTargs)>700 & Delay(1,AllTargs)<900, Delay(2,AllTargs)=3; end; catch; end;
                    SubID(AllTargs) = s;
                    Targ(AllTargs) = trigs(i);
                    
                end   
                
            end
            
                
        end


    
conds =[1:3; 4:6]; 
Data1 = [Delay(2,SubID==s).*(Targ(SubID==s)); Accuracy(SubID==s); RT(SubID==s)]; 
Data1(1,Data1(1,:)==8)=1; Data1(1,Data1(1,:)==16)=2; Data1(1,Data1(1,:)==24)=3; 
Data1(1,Data1(1,:)==9)=4; Data1(1,Data1(1,:)==18)=5;Data1(1,Data1(1,:)==27)=6; 
save([Datafolder num2str(subjID{s}) 'ModellingData'],'Data1');

end
RT(RT<-600)=NaN; 
save([Datafolder 'BehaviouralOutput'],'SubID', 'Accuracy', 'RT', 'Score', 'Block', 'Delay','Targ', 'Resp', 'SubTrials');

Data=[Delay(2,:).*(Targ); Accuracy; RT];
Data(1,Data(1,:)==8)=1; Data(1,Data(1,:)==16)=2; Data(1,Data(1,:)==24)=3; 
Data(1,Data(1,:)==9)=4; Data(1,Data(1,:)==18)=5;Data(1,Data(1,:)==27)=6; 
save([Datafolder 'ModellingData'],'Data');  
   
            