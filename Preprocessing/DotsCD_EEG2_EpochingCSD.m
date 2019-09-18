%% Delays Study - Ciara Devine
%% Pre Process Data - Average Reference, Baseline Correct, Artifact Rejection and Create final ERP matrices

clear all
close all



path = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\FilteredData\';
matfolder = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\PreProcessedData\';
Files = { '_NoHPF_NoLPF' '_NoHPF' '_HPF_0_05' '_HPF_0_15' '_HPF_0_25' '_HPF_0_35' '_HPF_0_45' '_HPF_0_55' '_HPF_0_65'}; % final part of file name changes depending on which HPF is being used.
SubFold = {'NoLPF_NoHPF\' 'LPF_NoHPF\' 'LPF_HPF_0_05\' 'LPF_HPF_0_15\' 'LPF_HPF_0_25\' 'LPF_HPF_0_35\' 'LPF_HPF_0_45\' 'LPF_HPF_0_55\' 'LPF_HPF_0_65\'}; % Which subfolder to look in depending on which high pass filter being used
Files2 = {  '_NoHPF_NoLPF' '_NoHPF' '_0_05HPF' '_0_15HPF' '_0_25HPF' '_0_35HPF' '_0_45HPF' '_0_55HPF' '_0_65HPF'}; % final part of file name changes depending on which HPF is being used.
dotsfolder = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\BDFs\';
badchansfilename='Badchans_dots';
allsubj=  { 'GD'  'KM'  'OM'  'AF' 'SB' 'DH' 'SM' 'SD' 'HM' 'KW' 'DM' 'AE' 'ST' 'KL' 'SY' 'WR' 'JE' 'CM' 'EC' 'DPM' 'TS' 'SCM'}; % In order of when they were tested
blocks = {[1:5] [1:6] [1:5] [1:10] [1:10] [1:10] [1:10] [1:10] [1:2 4:10] [1:9] [1:9] [1:8] [1:8] [1:7] [1:8] [1:9] [1:8] [1:8] [1:8] [1:8] [1:8] [1:8]}; % For HM Block 3 got split into block 3 & 4 due to recording issue

% %% Define bad channels per person
% try load([matfolder 'Badchans_dots']);
% catch
%     badchans={{[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]} {[]}};
% end
%  if length(badchans)<length(allsubj)
%      for ss=length(badchans):length(allsubj)
%      badchans{ss}={[]};
%      end
%  end
%      
%% LOAD BADCHANS
try
    load([matfolder badchansfilename]);
    CSD_on=1; 
    checkbadchans=0; 
catch
    for s=1:length(allsubj)
    badchans{s}=[];
    CSD_on=0; % If no pre existing bad chans file then don't do CSD on this run as not necessary for checking artefacts
    checkbadchans=1; % 
    end
end

%% Define epochs
fs=512; % sample rate
ts800 = round(-0.8*fs: (.8+2)*fs);   % time in samples: ERP epoch length (-0.5s - 2s) converted into samples. So ts = -256-1024 samples.
ts1200 = round(-0.8*fs: (1.2+2)*fs);
ts1600 = round(-0.8*fs: (1.6+2)*fs);
t800 = ts800*1000/fs; % t = time in ms. Time in ms is calculated using 'sample'
t1200 = ts1200*1000/fs;
t1600 = ts1600*1000/fs;
tt={t800 t1200 t1600};
tts={ts800 ts1200 ts1600};

slock_winds=tts; slock_wind=tt;
rlock_winds = round(-.6*fs:.6*fs); %
rlock_wind = rlock_winds*1000/fs;


%% Other important info

TargCodes = [8 9 12 13 16 17];
Delays = [1 1 2 2 3 3];

BLint = [600 700];   % baseline interval in ms

chanlocs = readlocs('cap128.loc'); % loads in the channel coordinates, necessary for plotting topographies and for channel interpolation
chanlocs = chanlocs(1:128);

VEOGchans = [131 132]; % Subject 13 is missing VEOG data... Blinks will be detected as artifacts in EEG
ARchans=1:128;
artifth = 100; % cut off/threshold for artifact rejections that will be applied to scalp channels. I.e. if an epoch exceeds 100uV then it is counted as an artifact.
artifact_th_v = 200; 

clear files mat
% checkbadchans = 0; % Set this to one if you want to just look at artifacts without CSDing and compiling output matrices etc etc! 
% CSD_on=1;
%% Generate G & H matrices for CSD transformation
SplineSettings = [2:7]; % this denotes the flexibility value in the CSD transformation (can be between 2-10)

if CSD_on
try load('C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\BDFs\CSD_GH');
catch
% First Generate M
for spline=1:SplineSettings(end)
if spline>1
E = textread('C:\Users\devineca\Documents\MATLAB\CSDtoolbox\chans.asc','%s');
M = ExtractMontage('C:\Users\devineca\Documents\MATLAB\CSDtoolbox\cap128.csd',E);  % reading in the montage for the CSD toolbox
MapMontage(M); close all;
% Now generate G & H
[g, h]=GetGH(M,spline);
G{spline}=g; H{spline}=h;
else
G{spline}=NaN; H{spline}=NaN;
end
end
save('C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\BDFs\CSD_GH', 'G', 'H')
end
end
% G1=G; H1=H;
% load('C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\DotsDelays\CSD_coords')

%% LOOP 1: HIGH PASS FILTER LOOP
for HPF=[3];

%% LOOP 2:SUBJECTS LOOP
for s=1:length(allsubj); 
        
        
%% --Pre allocate / blank matrices for each subject-- %%
numtr=0; sessid = []; blockid=[]; blockid_ws=[];allrespLR=[];allRT=[]; allRTsamp=[]; allTrig=[]; 
delay=[]; delay2=[]; Correct_Error=[]; TrialLen=[];Artifs=[]; Blinks=[];
erp = single(NaN(128,length(ts1600)+310, 600)); 
if CSD_on
    erp_csd = single(NaN(128,length(ts1600)+310, 600));
end
NumTrials(s)=0;
f = 0; % set file counter to 0 for each subject



%% LOOP 3: BLOCK LOOP

for b=blocks{s}
f = f+1;
Noep=[];


%% Load in Raw Data 

filename=[path SubFold{HPF} allsubj{s} num2str(b) Files2{HPF}];
load(filename);

dotsfile = [dotsfolder allsubj{s} 'Dots'];
clear par motionenergy_trials
try load(dotsfile); catch; end;

%% Resample if sampling rate not 512
if EEG.srate>512 || EEG.srate<512,
EEG = pop_resample(EEG, 512);
end



%% Separate VEOG & EEG
VEOGData=EEG.data(VEOGchans,:);
EEG.data=EEG.data(1:128,:);


%%  INTERPOLATE BAD CHANNELS
EEG.chanlocs = chanlocs;EEG.nbchan=128;
if ~isempty(badchans{s})
EEG=eeg_interp(EEG, badchans{s}{:},'spherical'); %
end

%%  AVERAGE-REFERENCE
numchan = 128;
% Applied to the whole continuous data (safe to do this now after interpolation - otherwise bad channels included in average reference):
EEG.data = EEG.data - repmat(mean(EEG.data([1:128],:),1),[numchan,1]);

%% Extract triggers & timings
clear trigs stimes RT artifact artifact_v 
for i=1:length(EEG.event)
trigs(i)=EEG.event(i).type;
stimes(i)=round(EEG.event(i).latency);
end

targtrigs = []; resptrigs = []; stimtrigs = []; Offset=[]; fixtrigs=[]; 
for n=1:length(trigs)-3
if trigs(n) == 24 && ismember(trigs(n+1), TargCodes) && ismember(trigs(n+2), [20 25]) && trigs(n+3)==5
targtrigs = [targtrigs n+1];
stimtrigs = [stimtrigs n];
resptrigs = [resptrigs n+2];
Offset = [Offset n+3];
fixtrigs = [fixtrigs n-1];
elseif trigs(n) == 24 && ismember(trigs(n+2), TargCodes) &&  ismember(trigs(n+1), [20 25]) && trigs(n+3)==5;  % Don't include a stim trigger if it is the final trigger in trigs ie not followed by any target
targtrigs= [targtrigs n+2];
stimtrigs = [stimtrigs n];
resptrigs = [resptrigs n+1];
Offset = [Offset n+3];
fixtrigs = [fixtrigs n-1];
elseif trigs(n) == 24 && ismember(trigs(n+1), TargCodes) && trigs(n+2)==5;  % Don't include a stim trigger if it is the final trigger in trigs ie not followed by any target
targtrigs= [targtrigs n+1];
stimtrigs = [stimtrigs n];
resptrigs = [resptrigs NaN];
Offset = [Offset n+2];
fixtrigs = [fixtrigs n-1];
end
end

stimtrigs = stimtrigs(1:length(targtrigs)); % only accept stimtrigs for which there is a corresponding targtrig

try
par.resptime{b}(par.resptime{b}(:)==0)=nan;
par.respLR{b}(par.respLR{b}(:)==0)=3;
catch
end


%% LOOP 4: TRIAL/EPOCH LOOOP
for n=1:length(targtrigs) 

    %% Retrieve timing of events for the current trial 
    stimtime = stimes(stimtrigs(n)); % stimtime = time of stimulus onset 
    evtime = stimes(targtrigs(n)); % evtime = time of target change
    Offtime = stimes(Offset(n)); % time stimulus leaves screen
    try fixtime = stimes(fixtrigs(n)); catch; end; % time of fixation dot
    try resptime = stimes(resptrigs(n)); catch; resptime=NaN; end
    Delay2 = evtime-stimtime;
    if Delay2<500, Delay = 1; elseif Delay2>600 && Delay2<700, Delay=2; elseif Delay2>800, Delay=3; end;
    TargLen = Offtime - evtime;
    FixLen = stimtime-fixtime;

    %% Create Epoch
    clear ep epBL; ep=[]; VEOGep=[]; EvOnsetSamp=[];
    clear t ts; t=tt{Delay}; ts=tts{Delay};

    if stimtime>abs(ts(1)) && Offtime+310<size(EEG.data,2) % make sure that the epoch doesn't go outside of the data

    ep = double(EEG.data(:,stimtime+ts(1):stimtime+ts(end) + 310));
    VEOGep = VEOGData(1,stimtime+ts(1):stimtime+ts(end)+310)-VEOGData(2,stimtime+ts(1):stimtime+ts(end)+310);
    epfix=double(EEG.data(:,fixtime-200:fixtime+400));
    EvOnsetSamp=find(ts==Delay2); % Find the sample point corresponding to 800ms post stimulus onset by subtracting 800 from all time points in t and getting the min of the abs values. Evidence Onset corresponds to the ~565th sample of the epoch ie 1100ms after the beginning of the epoch (epoch begins at -300).
    
    %% Calculate RT
    try
    thisRT = resptime - evtime;
    thisRTSamp = EvOnsetSamp+thisRT;
    valid_r=1;
    catch
    thisRT=NaN;
    thisRTSamp = NaN;
    valid_r=0;
    end
    end



    if thisRT<-(600/(1000/512)) || thisRT>(2000/(1000/512)) % Don't allow RTs through unless they are greater than -600ms & less than 2000ms (on some trials the trial may last longer than 2 seconds - this needs to be accounted for).
    ep=[];
    Noep(n)=1; 
    else
    Noep(n)=0;
    end


    %% Baseline Correction For Purpose of AR Rejection
    if ~isempty(ep)
    numtr=numtr+1;

    BLamp =  mean(ep(:,find(t>BLint(1) & t<BLint(2))),2);
    epBL = ep- repmat(BLamp,[1,size(ep,2)]);
    %eprBL = epr - repmat(BLamp,[1,size(epr,2)]);

%% ---------------------------------------------- CHECKING FOR ARTIFACTS STIMULUS LOCKED ------------------------------------------- %%

%% Create Windows for evaluating artifacts separately in the baseline period, delay period and main epoch perio


    if valid_r && thisRT>0;
    PreStimWind1=[]; PreStimWind1=find(ts==-400):find(ts==-200);
    PreStimWind2=[]; PreStimWind2=find(ts==-200):find(ts==0);
    EPWind=[]; EPWind = find(ts==0): EvOnsetSamp+thisRT+50;
    PostRTWind = []; PostRTWind = EvOnsetSamp+thisRT:length(epBL);
    else
    PreStimWind1=[]; PreStimWind1=find(ts==-400):find(ts==-200);
    PreStimWind2=[]; PreStimWind2=find(ts==-200):find(ts==0);
    EPWind=[]; EPWind = find(ts==-50): length(epBL);
    PostRTWind = []; PostRTWind = nan;
    end

    if EPWind(end)>length(epBL)
    EPWind = EPWind(EPWind<(length(epBL)+1));
    end
    %% Check & Record artifacts
    Artifs.PreStim1(:,numtr)= max(abs(epBL(ARchans,PreStimWind1)),[],2)>artifth;
    Blinks.PreStim1(numtr) = abs(max(VEOGep(PreStimWind1))-min(VEOGep(PreStimWind1)))>artifact_th_v;
    Artifs.PreStim2(:,numtr)= max(abs(epBL(ARchans,PreStimWind2)),[],2)>artifth;
    Blinks.PreStim2(numtr) = abs(max(VEOGep(PreStimWind2))-min(VEOGep(PreStimWind2)))>artifact_th_v;

    Artifs.Epoch(:,numtr)= max(abs(epBL(ARchans,EPWind)),[],2)>artifth;
    Blinks.Epoch(numtr) = abs(max(VEOGep(EPWind))-min(VEOGep(EPWind)))>artifact_th_v;


    if valid_r && thisRT>0
    Artifs.PostRTWind(:,numtr) = max(abs(epBL(ARchans,PostRTWind)),[],2)>artifth;
    Blinks.PostRTWind(numtr) = abs(max(VEOGep(PostRTWind))-min(VEOGep(PostRTWind)))>artifact_th_v;
    else
    Artifs.PostRTWind(:,numtr) = zeros(length(ARchans),1);
    Blinks.PostRTWind(numtr) = 0;
    end


    %% Compile Output Matrices
    if ~checkbadchans
    % stimulus locked epoch
    erp(:,1:length(ep),numtr)= ep;
    erpfix(:,:,numtr)=epfix;
    %CSD
    if CSD_on
    erp_csd(:,1:length(ep),numtr) = CSD(ep, G{4}, H{4});
    erpfix_csd(:,:,numtr) = CSD(epfix, G{4}, H{4});
    end


    allTrig(numtr) = trigs(targtrigs(n)); % find the target in trigs that corresponds to this epoch (n) Ie an
    try allrespLR(numtr)=trigs(resptrigs(n)); catch; allrespLR(numtr)=NaN; end
    allRT(numtr) = thisRT;
    allRTsamp(numtr) = thisRTSamp;
    blockid(numtr) = b;
    delay(numtr) = Delay;
    delay2(numtr) = Delay2;
    Targlen(numtr)=TargLen;
    Fixlen(numtr)=FixLen;

    if (ismember(allTrig(numtr), [8 12 16]) && allrespLR(numtr)==20) || (ismember(allTrig(numtr), [9 13 17]) && allrespLR(numtr)==25); Correct_Error(numtr)=1;
    elseif (ismember(allTrig(numtr), [8 12 16]) && allrespLR(numtr)==25) || (ismember(allTrig(numtr), [9 13 17]) && allrespLR(numtr)==20); Correct_Error(numtr)=0;
    else Correct_Error(numtr)=2; % Code misses
    end

    %% Match up Dots data with EEG & behavioural data from BDFS
    % For a match the delay, target type, response type and RT must be the same. There is a margin of
    % error in RTs recorded in the BDF files vs mat files of ~ 13-17ms 
    % Some trials for certain subjects won't have a match due to the appropriate mat file being
    % missing. When this happens dots data will be recorded as NaNs 

    %%% All Possible matches
    % 1. Perfect Match: Delay, target type, response type and response time all match 
    % 2. Imperfect Match 2(Genuine Miss trial): Delay & target type match but response=3 and response time=nan.  
    % 3. Imperfect Match 1 (Fast RTs or very slow RTs): Delay & target type matched but response & RT recorded only in BDF. Matfile equivalents=3 & nan.

    %%% All possible mismatches 
    % 1. Regular: delay &/or target type not matched.
    % 2. Regular2: delay & target & response type matched but RT mismatched.
    % 3. Regular3: delay & target matched but response type ~=3 in mat file & RT ~=nan (ie not a fast trial match).

    try


    % Convert to a common coding frame ie 1= (targcodes) 8 12 16 & respcodes(20
    % 2=(targcodes) 9 13 17 & respcodes (25)
    Targcodes2{1}=[8 12 16];
    Targcodes2{2}=[9 13 17];
    Respcodes=[20 25 5];
    if ismember(allTrig(numtr), Targcodes2{1}), matcode=1;
    else matcode=2;
    end

    % Change nans to 3 for consistency with mat file 
    if isnan(allrespLR(numtr)), thisresp=nan; 
    % code responses using 1 and 2 instead of 20/25
    else thisresp=find(Respcodes==allrespLR(numtr));
    end

    %%% Indices for Perfect Match trials 
    % Index trials where delay, target and response type match between mat and BDF 
    dotsidxA=[];dotsidxA = par.Delay{b}(:)==Delay & par.trialLR{b}(:)==matcode & par.respLR{b}(:)==thisresp;
    % Index trials for which there is a match for the current trial RT
    dotsidxB1=[];dotsidxB1 = ismember(round((round(par.resptime{b}(:)*1000) - round(allRT(numtr)*(1000/fs)))), [13:18]); 

    %%% Indices for potential Imperfect Match trials 
    % Index trials where delay & target match but
    % response type doesn't (some of these trials are genuine misses but some may have been trials where responses weren't recorded in the matfile eg ''clicked too soon'' 
    dotsidxA2=[];dotsidxA2 = par.Delay{b}(:)==Delay & par.trialLR{b}(:)==matcode & par.respLR{b}(:)==3;
    % Index trials where no responsetime is recorded. Some of these will be genuine misses but others may be 'clicked too soon' trials 
    dotsidxB2=[];dotsidxB2=isnan(par.resptime{b}(:)); 

    %Genuine misses can be identified where no response (ie nan) is recorded in the BDF file data


    if allRT(numtr)<75
        stophere=1;
    end
    %%% First look for trials in which RT was faster than 150ms according to BDF but recorded as a
    % nan in the matfile. Hone in on those trials that are already matched for delay and target type

    %%%% First check that the current trial is a valid response (ie not a miss)    
    if ~isnan(allRT(numtr)) & ~isnan(allrespLR(numtr)) 
        % If there is a valid response, check for trials in the matfile that match the present trial.
        % If there are no matches then record nans. %%%(this may occur sometimes due to issues with synchronising mat files and BDF files at time of recording)
        if sum(dotsidxA & dotsidxB1)==0 & sum(dotsidxA2 & dotsidxB2)==0
        TrialMatch(numtr,1)=nan;  
        TrialMatch(numtr,2)=nan;
        TrialMatch(numtr,3)=b;
        Matches{numtr}=nan;
        Dots{numtr}=nan;

        % If there is a match (either a perfct match or an imperfect match (described above)) 
        elseif sum([(dotsidxA) & (dotsidxB1)])~=0 | sum([(dotsidxA2) & (dotsidxB2)])~=0 
            Matches{numtr}=find([dotsidxA & dotsidxB1] | [dotsidxA2 & dotsidxB2]); % Get indices for all possible matches 
            [srt srtidx]=sort(abs(Matches{numtr}-n));
            TrialMatch(numtr,1)=Matches{numtr}(srtidx(1)); % Choose the match that is closest in ordering/sequence to the BDF trial
            TrialMatch(numtr,2)=TrialMatch(numtr,1)+sum(par.numtrials(1:b-1));
            TrialMatch(numtr,3)=b;
            Dots{numtr}=motionenergy_trials{TrialMatch(numtr,1)+sum(par.numtrials(1:b-1))};
            clear dotdirection; dotdirection=cumsum(motionenergy_trials{TrialMatch(numtr,1)+sum(par.numtrials(1:b-1))});
            if dotdirection<0 & matcode==2
                stophere=1;
            elseif dotdirection>0 & matcode==1
                stophere=1;
            end
        end

    %                             
    %%% If there was no valid response ie participant incurred a ''miss''
    elseif isnan(allRT(numtr)) & isnan(allrespLR(numtr)) % Look also for trials where there genuinely was no response and hence only delay and target type will be matched  
        % Check for trials in matfile that provide an imperfect match with the current trial
        % If no imperfect match, record nans.
        if sum(dotsidxA2 & (dotsidxB2))==0
        TrialMatch(numtr,1)=nan; 
        TrialMatch(numtr,2)=nan;
        TrialMatch(numtr,3)=b;
        Matches{numtr}=nan;
        Dots{numtr}=nan;

        % if imperfect matches are found
        elseif sum(dotsidxA2 & dotsidxB2)~=0
            Matches{numtr}=find(dotsidxA2 & dotsidxB2);
            [srt srtidx]=sort(abs(Matches{numtr}-n));
            TrialMatch(numtr,1)=Matches{numtr}(srtidx(1));
            TrialMatch(numtr,2)=TrialMatch(numtr,1)+sum(par.numtrials(1:b-1));
            TrialMatch(numtr,3)=b;
            Dots{numtr}=motionenergy_trials{TrialMatch(numtr,1)+sum(par.numtrials(1:b-1))};
            Dots{numtr}=motionenergy_trials{TrialMatch(numtr,1)+sum(par.numtrials(1:b-1))};
            %clear dotdirection; dotdirection=cumsum(motionenergy_trials{TrialMatch(numtr,1)+sum(par.numtrials(1:b-1))});
            % Debugging: Check that the direction of dots is consistent with the matched trial 
            if dotdirection<0 & matcode==2
                stophere=1;
            elseif dotdirection>0 & matcode==1
                stophere=1;
            end
        end
    end

    %                             %%%%%%%% Code for troubleshooting - 
    %                             stop at trials that have more than 1 match and ensure
    %                             that code is choosing the appropriate match
    %                             if length(Matches{numtr})>1 || length(Matches{numtr})<1
    %                                 stophere=1;
    %                             end

    %                             % Stop if any of the following irregularities occur 
    %                             if numtr>1 & n>1
    %                                 if ~ismember((TrialMatch(numtr,1)-TrialMatch(numtr-1,1)), [1 1-par.numtrials(b)]) ...
    %                                     & ~Noep(n-1) & ~(length(targtrigs)<par.numtrials(b));
    %                                     stophere=1;
    %                                 end
    %                             end
    %%%%%%%%%%%%%%%%%%%%%%%%%

    catch  % some subjects did not have a full set of mat trials (eg KL) - allocate nans instead 
    TrialMatch(numtr,1)=nan;
    TrialMatch(numtr,2)=nan;
    TrialMatch(numtr,3)=b;
    Matches{numtr}=nan;
    Dots{numtr}=nan;

    end
    end

    end
end
end


%% Choose channels to add to the list for interpolation in future
AnyArt=[]; AnyBlink=[]; NewArChans=[]; ArtNoBlink=[];

AnyArt=(Artifs.Epoch)>0;
AnyBlink= (Blinks.Epoch)>0;
ArtNoBlink = sum(AnyArt(:,~AnyBlink),2);


if s==11
NewArtChans=find(ArtNoBlink>60)';
else
NewArtChans=find(ArtNoBlink>27)';
end

badchans{s}{:} = [badchans{s}{:} NewArtChans];

numtrials = sum((sum(AnyArt,1) + AnyBlink)<1);

%%  Plots histogram to show frequency of artifacts in each channel for each subject/session
figure; hold on;
subplot(1,3,1);
bar(1:128, sum(Artifs.PreStim2(:, ~AnyBlink),2)); hold on; box on;
title([allsubj{s} 'BL NumTrials = ' num2str(numtrials) ' Blinks=' num2str(sum(AnyBlink))]);

subplot(1,3,2);
bar(1:128, sum(Artifs.Epoch(:, ~AnyBlink),2)); hold on; box on;
title([allsubj{s} 'StimOn: 100ms PostRT']);

subplot(1,3,3);
bar(1:128, sum(Artifs.PostRTWind(:, ~AnyBlink),2)); hold on; box on;
title([allsubj{s} 'Post RT']);

savefig([matfolder SubFold{HPF} allsubj{s} 'Artifacts']);

close all;

%% Trim Away Nans
if ~checkbadchans
erp=single(erp(:,:,1:numtr));   erpfix=single(erpfix); 
if CSD_on
erp_csd=single(erp_csd(:,:,1:numtr)); erpfix_csd=single(erpfix_csd);
end

allRT=single(allRT); allrespLR=single(allrespLR); delay=single(delay);
delay2=single(delay2); sessid=single(sessid); blockid=single(blockid);
%TrialLen=single(TrialLen);TargLen=single(TargLen);
allRTsamp = single(allRTsamp); 

%% Save Data
save([matfolder SubFold{HPF} allsubj{s} 'DelaysERP_NoCSD' Files{HPF}],'erp','erpfix', 'tt', 'tts', 'rlock_winds', 'rlock_wind', 'Artifs','Blinks', 'fs', 'allRT', 'allRTsamp','allTrig', 'allrespLR','blockid', 'delay', 'delay2','Correct_Error') %,'erp_hpf', 'erp_detrend','erp_hpfdetrend',

if CSD_on
clear erp erpfix
erp=erp_csd; erpfix=erpfix_csd; 
save([matfolder SubFold{HPF} allsubj{s} 'DelaysERP_CSDSpline4_' Files{HPF}],'erp', 'erpfix','tt', 'tts', 'rlock_winds',  'rlock_wind', 'slock_winds', 'slock_wind','Artifs', 'Blinks','fs', 'allRT', 'allRTsamp','allTrig', 'allrespLR','blockid','delay', 'delay2','Correct_Error') %,'erp_hpf', 'erp_detrend','erp_hpfdetrend',
end
end

save([matfolder 'Badchans_dots'], 'badchans');
save([matfolder SubFold{HPF} allsubj{s} 'Dots'], 'Dots', 'TrialMatch', 'Matches');
clear Dots TrialMatch Matches;
end

end


%%%%CSD TRANSFORMATION%%%% e.g. Keysers et al. 2015

%Critical problems in EEG: EEG has poor spatial resolution and electrical signals recorded at the scalp generally represent blurred and distorted copies
%of the original cortical activity. Due to volume conduction, eletrical activity from dipoles in the brain spreads out across the scalp and is picked up by electrodes located across the scalp as opposed to just directly above the sites of cortical activity. With high levels of resistance
% in the skull, electrical activity at the scalp level tends to be highly distributed and the head/skull essentially acts as a low pass spatial filtern transmitting broad as opposed to focal spatial patterns of activity to the scalp. As a result EEG signals are generally composed of activity from a mixture of cortical sources.
% Electrical activity at one location therefore can be contaminated by activity from sources in the brain that are highly distributed and located far from the location of interest. In other words, activity recorded at each electrode does not represent specific activity located
% from local sources in the brain (ie directly beneath the electrode) but instead represents volume conducted activity coming froms spatially distributed sources.

%Another limitation in EEG research is that EEG activity is reference dependent (ie electrical potentials require some point of reference) and therefore scalp based EEG measures cannot be understood in isolation from their reference. The choice of reference scheme affects the resulting
%waveforms in EEG analyses. Thus the results of EEG analyses depend on the reference used. ERP waveforms therefore always reflect contributions both from the active and the reference sites and waveforms can look different with different references.

%Current Source Density (CSD, aka, Surface Laplacian (SL) or Scalp Current Density (SCD)) is a method of signal enhancement  or ''spatial deblurring'' designed to improve the frequency resolution of EEG data and counteract the problems outlined above.
%Effectively CSD filters out spatially broad features of the EEG data in order to counteract the distortions caused by volume conduction, reduce the contribution of distant sources of cortical activity at each eletrode site and to facilitate the spatial separation of components.
% It is a spatial band pass filter that hones in on high spatial frequency activity and attenuates low spatial frequency activity thereby preserving activity which is visible in small clusters of neurons but not activity that affects/spreads to larger areas of the scalp.
% It thereby minimises volume conduction distortions and improves/sharpens topographical localisation making CSD data more spatially distinct.

%The method involves converting voltage data into current density data which, rather than measuring the potential for charges to flow between two points (ie voltage), measures/estimates the flow of current radially towards and out of the scalp from underlying neural tissue at each electrode site/scalp location (ie current).
% (Voltage is electrical potential (or difference in charge) between two points in a circuit and represents the force that movement electrons through a circuit.
% Current on the other hand refers to the actual flow of the electronsthrough the circuit. It is measured in terms of the number of electrons flowing through a point in the circuit over a unit of time.
% Current density in turn refers to the amount of current per unit of area).

%Computationally, CSD works by calculating the ''second spatial derivative'' of the voltage distribution over the scalp. In simpler terms it involves summing the spatially weighted potential gradients directed at/contributing to each electrode site.

% LOCAL METHOD: At its most basic level of implementation, the procedure involves computing the difference between the potential at each electrode site and the average potential at the sites of its nearest neighbours. Eg the 'local' Hjorth method (Hjorth, 1975)
% Neighbouring sites are weighted by their inverse distance from the potential site (using a differentiation grid that is independent of the actual data). Computing this difference is called calculating the first derivative. Repeating this process again using the first derivative values gives you the second derivative values.
% The problem with a local method is that ''nearest neighbours'' are not as plentiful for peripheral electrodes.

% GLOBAL METHOD: Global approaches compute CSD using interpolation functions. Spherical spline interpolation is a popular method (Perrin et al 1989).
% Interpolation is used so as to gauge a continuous measure or estimate of the voltage distribution over the scalp at each moment in time. Interpolation uses data (voltage measurements) at the actual eletrode sites in order to estimate and provide a continuous projection of the voltage at missing/intermediate sites on the scalp
% The interpolated voltage distribution values are then used as the basis for calculating the second spatial derivative values (which become the new CSD-transformed EEG data).
% More accurate interpolation can be carried out with higher density EEG because it is the measured data that is used in interpolation to provide a continuous projection of missing data on the scalp surface (ie  locations between eletrodes).

% Using this method based on spherical spline interpolation the outcome is that a weighted sum of activity at all eletrodes is subtracted from each eletrode.
% G and H matrices contain the weights that are applied to data (the first derivative is calculated using G and second derivative using H). These matrices/weightings are based on interelectrode, spherical distances
% The scalp is modelled as a sphere and distance between electrodes is measured as the cosine of the angle between surface locations. Unlike nearest neighbour approaches, all recording sites affect the laplacian but their influence dependings on their location and interelectrode distances.

%Another advantage of using current density instead of voltage is that current density is naturally defined at single points and therefore does not require a reference. In other words, CSD data are reference free estimates of radial current flow.

%A note on spline flexibility: Spline flexibility = degree to which the spherical spline function used to fit the observed data can be bent to best fit the actual data. Small m value = more flexible. M=2 is the most flexible.
%With higher flexibility eg m =2 or m=3 surface laplacian estimates are sharper and more enhanced whereas with lower flexibility in the splines eg m=4 or m=5 then laplacian estimates are smoother and more gradual
% Which spline setting is better?? To change spline setting type ''open GetGH'' in command window


% 
%  %
%  nn=n;
%  
%  for check=[1 2 4:50] % 3 options - check match for the current trial and failing that check for a match with the previous or subsequent trial
%      try
%          if TrialMatch(numtr)==0 && ismember(check,[2 4:50])  % Allow for the possibility that matched trials have become disordered
%              n=nn+(check-3); % check-3 =-1 (check previous trial), check-3=1 (check subseqential trial)
%          else
%          end
%      catch
%      end
%      if n==0, % Hack to prevent n ever being zero ie for the first trial in a block
%          n=1; end
%      
%      % Check that the trials in the mat file and those
%      % in the BDF match in terms of delay, target type
%      % and response type
%      if (Delay==par.Delay{b}(n)) && ismember(allTrig(numtr),Targcodes2{par.trialLR{b}(n)})...
%              &&  (ismember(allrespLR(numtr),Respcodes(par.respLR{b}(n))) || (isnan(allrespLR(numtr)) ...
%              && par.respLR{b}(n)==3) || (ismember(allrespLR(numtr), [20 25]) && par.respLR{b}(n)==3));
%          % Now check that the RTs match
%          if abs(round(allRT(numtr)*(1000/fs)) - (par.resptime{b}(n)*1000))<20
%              TrialMatch(numtr)=1;
%              dotidx = (par.numtrials(b)*(b-1))+n;
%              Dots{numtr}=motionenergy_trials{dotidx};
%              
%          elseif round(allRT(numtr)*(1000/fs))<152 && isnan(par.resptime{b}(n)) % RT was recorded in structure 'par' as 0 if faster than 150ms & then changed to nan. We are still interested in these fast responses though!!
%              TrialMatch(numtr)=1;
%              dotidx = (par.numtrials(b)*(b-1))+n;
%              Dots{numtr}=motionenergy_trials{dotidx};
%              
%          elseif isnan(allRT(numtr)) && isnan((par.resptime{b}(n)*1000)) % Trials without responses can still be a match
%              TrialMatch(numtr)=1;
%              dotidx = (par.numtrials(b)*(b-1))+n;
%              Dots{numtr}=motionenergy_trials{dotidx};
%          else
%              % If there is no match, check that the number
%              % of trials in the mat and BDF files match
%              % If there are fewer trials in the BDF than the
%              % mat then check if the current BDF trial
%              % matches the next mat trial
%              TrialMatch(numtr)=0;
%              Dots{numtr}=nan;
%          end
%      else
%          TrialMatch(numtr)=0;
%          Dots{numtr}=nan;
%      end
%  end
%  
%  catch
%      TrialMatch(numtr)=0;
%      Dots{numtr}=nan;
%      end
%      





