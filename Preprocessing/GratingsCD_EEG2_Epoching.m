%% Delays Study - Ciara Devine
%% Pre Process Data - Average Reference, Baseline Correct, Artifact Rejection and Create final ERP matrices

clear all
close all



path = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\FilteredData\'];
matfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\EpochedData\'];
Files = {'_NoHPF' '_HPF_0_05' '_HPF_0_15' '_HPF_0_25' '_HPF_0_35' '_HPF_0_45' '_HPF_0_55' '_HPF_0_65'}; % final part of file name changes depending on which HPF is being used.
SubFold = {'LPF_NoHPF\' 'LPF_HPF_0_05\' 'LPF_HPF_0_15\' 'LPF_HPF_0_25\' 'LPF_HPF_0_35\' 'LPF_HPF_0_45\' 'LPF_HPF_0_55\' 'LPF_HPF_0_65\'}; % Which subfolder to look in depending on which high pass filter being used
Files2 = {'_NoHPF' '_0_05HPF' '_0_15HPF' '_0_25HPF' '_0_35HPF' '_0_45HPF' '_0_55HPF' '_0_65HPF'}; % final part of file name changes depending on which HPF is being used.
badchansfilename=['BadChans_Gratings'];
flatchansfilename=['FlatChans_Gratings'];
allsubj =  { 'CJG'  'CR'  'AR' 'KD'  'TB'  'AOB'  'PM'  'CMG'  'SR'  'AB' 'GK' 'ED' 'SC' 'CB' 'CE' 'SH' 'RMC' 'JH' 'JR' 'ROC' 'KM' 'ROS'}; % In order of when they were tested; KD and JR: too many trials lost due to blinks - excluded!
sessions = {[1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:4] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5]};
blocks =  { {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:29; 30:39; 40:49}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50},{1:10; 11:20; 21:30; 31:40; 41:50}, {1:10; 11:20; 21:30; 31:40; 41:50} ,{1:10; 11:20; 21:30; 31:40; 41:50}};


%% LOAD BADCHANS
try
    load([matfolder badchansfilename]);
    CSD_on=1; 
    checkbadchans=0; 
catch
    load([matfolder flatchansfilename]); %Load file that contains list of flat channels
    CSD_on=0; % If no pre existing bad chans file then don't do CSD on this run - Check for artifacts first instead
    checkbadchans=1; % 
end

%% Define epochs

fs=512; % sample rate

ts800 = round(-0.4*fs: (.8+2.05)*fs);   % time in samples: ERP epoch length (-0.5s - 2s) converted into samples. So ts = -256-1024 samples.
ts1200 = round(-0.4*fs: (1.2+2.05)*fs);
ts1600 = round(-0.4*fs: (1.6+2.05)*fs);
t800 = ts800*1000/fs; % t = time in ms. Time in ms is calculated using 'sample'
t1200 = ts1200*1000/fs;
t1600 = ts1600*1000/fs;

tt={t800 t1200 t1600};
tts={ts800 ts1200 ts1600};



%% Other important info

TargCodes = [8 9];
RespCodes= [12 13];
Delays = [1 2 3]; % 800 1200 1600

BLint = [550 650];

chanlocs = readlocs('cap128.loc'); % loads in the channel coordinates, necessary for plotting topographies and for channel interpolation
chanlocs = chanlocs(1:128);

VEOGchans = [131 132];
ARchans=1:128;
artifth = 100; % cut off/threshold for artifact rejections that will be applied to scalp channels. I.e. if an epoch exceeds 100uV then it is counted as an artifact.
artifth2 = 200; % artifact threshold for period after response
artifact_th_v = 200; % cut off for artifact rejection that will be applied specifically to eye channels

checkbadchans = 0;
CSD_on=1;
%% Generate G & H matrices for CSD transformation
SplineSettings = [1:5]; % this denotes the flexibility value in the CSD transformation (can be between 2-10)

try load('C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\GratingsCD_Data\CSD_GH');
catch
    % First Generate M
    for spline=1:SplineSettings(end)
        if spline>2
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
    save('C:\Users\devineca\Desktop\Learning_data\BDFs\CSD_GH', 'G', 'H')
end
%G1=G; H1=H;

%load('C:\Users\devineca\Desktop\Learning_data\BDFs\CSD_coords');

%% Loop 1: High pass filters
for HPF=2
    %% Loop 2: Loop Through All Subjects
    for s= 1:length(allsubj); % loop through all subjects
        
        
        %% --------------- Pre allocate / blank matrices for each subject----------------- %%
        clear erp* Exclude
        numtr=0; sessid = []; blockid=[]; blockid_ws=[];allrespLR=[];allRT=[]; allTrig=[]; delay=[]; delay2=[]; Correct_Error=[];
        Artifs=[]; Blinks=[]; AnyArt=[]; AnyBlink=[]; allRTsamp=[];
        
        erp = NaN(128,length(ts1600), 2500, 'single');
        if CSD_on
            erp_csd = NaN(128,length(ts1600), 2500,'single');
        end
        NumTrials(s)=0;
        
        %% Loop 3: Loop through Sessions
        for d=1:length(sessions{s})
            
            
            %% Loop 4: Loop through all blocks
            f = 0; % set file counter to 0
            for b=blocks{s}{d}
                f = f+1;
                
                
                
                %% Load in Raw Data
                
                filename=[path SubFold{HPF} allsubj{s} num2str(b) Files2{HPF}];
                load(filename);
                
                
                %% Resample if sampling rate not 512
                if EEG.srate>512 || EEG.srate<512,
                    EEG = pop_resample(EEG, 512);
                end
                
                %% Separate VEOG & EEG
                VEOGData=single(EEG.data(VEOGchans,:));
                EEG.data=single(EEG.data(1:128,:));
                
                
                %% INTERPOLATE BAD CHANNELS
                EEG.chanlocs = chanlocs;EEG.nbchan=128;
                if ~isempty(badchans{s}{d}{:})
                    EEG=eeg_interp(EEG, badchans{s}{d}{:},'spherical'); %
                    NumBadChans(d)=length(badchans{s}{d}{:});
                end
                
                
                %% AVERAGE-REFERENCE
                numchan = 128;
                % Applied to the whole continuous data (safe to do this now after interpolation - otherwise bad channels included in average reference):
                EEG.data = EEG.data - repmat(mean(EEG.data([1:128],:),1),[numchan,1]);
                
                %% EXTRACT TRIGGERS AND TIMES
                clear trigs stimes RT artifact artifact_v %
                for i=1:length(EEG.event)
                    trigs(i)=EEG.event(i).type;
                    stimes(i)=round(EEG.event(i).latency);
                end
                
                
                targtrigs = []; resptrigs = []; stimtrigs = []; Offset=[]; fixtrigs=[]; Early=[];
                for n=1:length(trigs)-3
                    if trigs(n) == 24 && ismember(trigs(n+1), TargCodes) && ismember(trigs(n+2), RespCodes) && trigs(n+3)==5
                        targtrigs = [targtrigs n+1];
                        stimtrigs = [stimtrigs n];
                        resptrigs = [resptrigs n+2];
                        fixtrigs = [fixtrigs n-1];
                        Offset = [Offset n+3];
                    elseif trigs(n) == 24 && ismember(trigs(n+2), TargCodes) &&  ismember(trigs(n+1), RespCodes) && trigs(n+3)==5;  % Don't include a stim trigger if it is the final trigger in trigs ie not followed by any target
                        targtrigs= [targtrigs n+2];
                        stimtrigs = [stimtrigs n];
                        fixtrigs = [fixtrigs n-1];
                        resptrigs = [resptrigs n+1];
                        Offset = [Offset n+3];
                    elseif trigs(n) == 24 && ismember(trigs(n+1), TargCodes) && trigs(n+2)==5;  % Don't include a stim trigger if it is the final trigger in trigs ie not followed by any target
                        targtrigs= [targtrigs n+1];
                        stimtrigs = [stimtrigs n];
                        resptrigs = [resptrigs NaN];
                        fixtrigs = [fixtrigs n-1];
                        Offset = [Offset n+2];
                    end
                end
                
                stimtrigs = stimtrigs(1:length(targtrigs)); % only accept stimtrigs for which there is a corresponding targtrig
                
                
                %% LOOP 4: TRIAL/EPOCH LOOP
                
                for n=1:length(targtrigs) % n is looping through targets - or further down looping through epochs
                    
                    %% Retrieve timing of events for the current trial
                    evtime = stimes(targtrigs(n)); % evtime = time of target onset
                    stimtime = stimes(stimtrigs(n)); % stimtime = time of stimulus onset (ie before target).
                    offtime = stimes(Offset(n));
                    try fixtime = stimes(fixtrigs(n)); catch; end;
                    try resptime = stimes(resptrigs(n)); catch; resptime=NaN; end
                    Delay2 = evtime-stimtime;
                    if Delay2<500, Delay = 1; elseif Delay2>600 && Delay2<700, Delay=2; elseif Delay2>800, Delay=3; end;
                    TargTime = offtime - evtime;
                    try ITITime = stimtime-stimes(Offset(n-1)); catch; ITITime=NaN; end
                    FixTime = stimtime-fixtime;
                    
                    
                    %% Create Epoch
                    clear ep ep2 epBL epBL2; ep=[]; VEOGep=[]; EvOnsetSamp=[];
                    clear t ts; t=tt{Delay}; ts=tts{Delay};
                    
                    if stimtime>abs(ts(1)) && offtime<size(EEG.data,2)
                        ep = EEG.data(:,stimtime+ts(1):stimtime+Delay2+1024);
                        VEOGep = VEOGData(1,stimtime+ts(1):stimtime+Delay2+1024)-VEOGData(2,stimtime+ts(1):stimtime+Delay2+1024);
                        %epfix=double(EEG.data(:,fixtime-200:fixtime+400));
                        EvOnsetSamp=find(ts==Delay2); % Find the sample point corresponding to 800ms post stimulus onset by subtracting 800 from all time points in t and getting the min of the abs values. Evidence Onset corresponds to the ~565th sample of the epoch ie 1100ms after the beginning of the epoch (epoch begins at -300).
                        %% Calculate RT
                        try
                            thisRT = resptime - evtime;
                            thisRTSamp = EvOnsetSamp+thisRT;
                            valid_r=1;
                        catch
                            thisRT=NaN;
                            thisRTSamp=nan;
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
                        
                        
                        
                        %% --------------------- CHECKING FOR ARTIFACTS STIMULUS LOCKED --------------------- %%
                        
                        %% Create Windows for evaluating artifacts separately in the baseline period, delay period and main epoch period
                        
                        if valid_r && thisRT>0;
                            PreStimWind1=[]; PreStimWind1=find(ts==-200):find(ts==-100);
                            PreStimWind2=[]; PreStimWind2=find(ts==-100):find(ts==0);
                            EPWind=[]; EPWind = find(ts==0): EvOnsetSamp+thisRT+50;
                            PostRTWind = []; PostRTWind = EvOnsetSamp+thisRT:length(epBL);
                        else
                            PreStimWind1=[]; PreStimWind1=find(ts==-200):find(ts==-100);
                            PreStimWind2=[]; PreStimWind2=find(ts==-100):find(ts==0);
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
                            %erpfix(:,:,numtr)=epfix;
                            %CSD
                            if CSD_on
                                erp_csd(:,1:length(ep),numtr) = CSD(ep, G{4}, H{4});
                                %erpfix_csd(:,:,numtr) = CSD(epfix, G{4}, H{4});
                            end
                            
                            allTrig(numtr) = trigs(targtrigs(n)); % find the target in trigs that corresponds to this epoch (n) Ie an
                            try allrespLR(numtr)=trigs(resptrigs(n)); catch; allrespLR(numtr)=NaN; end
                            allRT(numtr) = thisRT;
                            allRTsamp(numtr) = thisRTSamp;
                            blockid(numtr) = b;
                            sessid(numtr)=d;
                            delay(numtr) = Delay;
                            delay2(numtr) = Delay2;
                            TargLen(numtr)=TargTime;
                            FixLen(numtr) = FixTime;
                            
                            if (ismember(allTrig(numtr), 8) && allrespLR(numtr)==12) || (ismember(allTrig(numtr), 9) && allrespLR(numtr)==13); Correct_Error(numtr)=1;
                            elseif (ismember(allTrig(numtr), 8) && allrespLR(numtr)==13) || (ismember(allTrig(numtr), 9) && allrespLR(numtr)==12); Correct_Error(numtr)=0;
                            else Correct_Error(numtr)=2; % Code misses
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
            
            if (s==4 && d==5) || (s==6 && d==3)
                NewArtChans=find(ArtNoBlink>49)';
            else
                NewArtChans=find(ArtNoBlink>27)';
            end
            
            badchans{s}{d}{:} = [badchans{s}{d}{:} NewArtChans];
            
            numtrials = sum((sum(AnyArt,1) + AnyBlink)<1);
            
            %%  Plots histogram to show frequency of artifacts in each channel for each subject/session
            figure; hold on;
            subplot(1,3,1);
            bar(1:128, sum(Artifs.PreStim2(:,sessid==d & ~AnyBlink),2)); hold on; box on;
            title([allsubj{s} num2str(d) 'BL NumTrials = ' num2str(numtrials) ' Blinks=' num2str(sum(AnyBlink))]);
            
            subplot(1,3,2);
            bar(1:128, sum(Artifs.Epoch(:,sessid==d & ~AnyBlink),2)); hold on; box on;
            title([allsubj{s} num2str(d) 'Stim On : 100ms PostRT']);
            
            subplot(1,3,3);
            bar(1:128, sum(Artifs.PostRTWind(:,sessid==d & ~AnyBlink),2)); hold on; box on;
            title([allsubj{s} num2str(d) 'Post RT']);
            
            savefig([matfolder SubFold{HPF} allsubj{s} num2str(d) 'Artifacts']);
            
            if checkbadchans,
                save([matfolder 'BadChans_GratingsNew'], 'badchans', 'NumBadChans');
            end
        end
        
        close all;
        
        %% Trim Away Nans
        if ~checkbadchans
            erp=single(erp(:,:,1:numtr));  % erpfix=single(erpfix);
            if CSD_on
                erp_csd=single(erp_csd(:,:,1:numtr)); %erpfix_csd=single(erpfix_csd);
            end
%             
            allRT=single(allRT); allrespLR=single(allrespLR); delay=single(delay);
            delay2=single(delay2); sessid=single(sessid); blockid=single(blockid);
            allRTsamp = single(allRTsamp);
            
            %% Save Data
            save([matfolder SubFold{HPF} allsubj{s} 'GratingsDelaysERP_NoCSD' Files{HPF}],'erp', 'tt', 'tts', 'Artifs', 'Blinks', 'fs', 'sessid','allRT', 'allRTsamp','allTrig','allrespLR','blockid', 'delay', 'delay2','Correct_Error', 'sessid') %,'erp_hpf', 'erp_detrend','erp_hpfdetrend',
            
            if CSD_on
                clear erp %erpfix
                erp=erp_csd; %erpfix=erpfix_csd;
                clear erpBL_csd4 erpBL2_csd4
                save([matfolder SubFold{HPF} allsubj{s} 'GratingsDelaysERP_CSDSpl4_' Files{HPF}],'erp', 'tt', 'tts', 'Artifs', 'Blinks', 'fs','sessid', 'allRT', 'allRTsamp','allTrig','allrespLR','blockid', 'delay', 'delay2','Correct_Error','sessid') %,'erp_hpf', 'erp_detrend','erp_hpfdetrend',
                
            end
        end
        
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












