clear all
close all

%% File and folder names

matfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\PreprocessedData\']; % Where to put outgoing data
plotfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\DotsPlots\'];
SubFold = {'LPF_NoHPF\' 'LPF_HPF_0_05\' 'LPF_HPF_0_15\' 'LPF_HPF_0_25\' 'LPF_HPF_0_35\' 'LPF_HPF_0_45\'  'LPF_HPF_0_55\'  'LPF_HPF_0_65\'}; % Which subfolder to look in depending on which high pass filter being used

matfolder2 = ['C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\ERPMatrices\']; % Where to put outgoing data
Files = {'_NoHPF' '_HPF_0_05' '_HPF_0_15' '_HPF_0_25' '_HPF_0_35' '_HPF_0_45' '_HPF_0_55' '_HPF_0_65'}; % final part of file name changes depending on which HPF is being used.

DotsFolder = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\BDFs\';
%%


allsubj=  { 'GD'  'KM'  'OM'  'AF' 'SB' 'DH' 'SM' 'SD' 'HM' 'KW' 'DM' 'AE' 'ST' 'KL' 'SY' 'WR' 'JE' 'CM' 'EC' 'DPM' 'TS' 'SCM'}; % In order of when they were tested
blocks = {[1:5] [1:6] [1:5] [1:10] [1:10] [1:10] [1:10] [1:10] [1:2 4:10] [1:9] [1:9] [1:8] [1:8] [1:7] [1:8] [1:9] [1:8] [1:8] [1:8] [1:8] [1:8] [1:8]}; % For HM Block 3 got split into block 3 & 4 due to recording issue


RespCode=[20 25];

fs = 512;
[B,A]=butter(4,6*2/fs); % creating filter for ERP waveforms
delays=[800 1200 1600];

HPF=2;
spline=4; 


%% Dots parameters
Fs = 100;

t_dots = {[10:10:2800],[10:10:3200],[10:10:3600]};
ts_dots={};
for delay = 1:3
    ts_dots{delay} = t_dots{delay}*Fs/1000;
end

ts_cohmo = [-0.750*Fs:0.750*Fs];
% ts_cohmo = [-1.000*Fs:0.750*Fs];
t_cohmo = ts_cohmo*1000/Fs;
trs = [-0.750*Fs:0.125*Fs];
tr = trs*1000/Fs;
            
%% Loop through subjects
SubID=[];Delay=[]; nn=0; Block=[]; RT=[]; Acc=[]; Targ=[]; Resp=[];
ERP=[]; ERPr=[]; FFTMotor=[]; FFTMotorR=[]; FFTOcc=[]; FFTOccR=[];FFTFrontal=[]; FFTFrontalR=[];
for s=1:length(allsubj)
    
    % Load in Data
    
   
    if ismember(spline, 2:7)
        load([matfolder, SubFold{HPF}, allsubj{s}, 'DelaysERP_CSDSpline',num2str(spline), '_', Files{HPF}]);
    elseif spline==1
        load([matfolder, SubFold{HPF}, allsubj{s}, 'DelaysERP_NoCSD', Files{HPF}]);
    end
    if s==16
    end
    filename = [matfolder SubFold{HPF}, allsubj{s} 'Dots'];
    try 
        load(filename); 
    catch 
    end;
        allRTsamp=allRTsamp(1:size(erp,3));
    slock_winds = tts;
    slock_wind = tt;
    
    rlock_winds=-1*fs:.6*fs;
    rlock_wind=rlock_winds*(1000/fs);
    erp=erp(:,1:length(slock_winds{3}),:);
    erp(erp==0)=NaN;    
    
   %% Exclude artifact trials
    
    AnyBlinks = []; AnyArts=[]; Exclude=[]; ZERPEv=[]; ZERPRT=[];
    
    AnyArts = (Artifs.Epoch)>0;
    Interp=find(sum(AnyArts,2)>30); % Check for really bad single channels that should have been interpolated
    erp(Interp,:,:)=nan; %erpr(Interp,:,:)=nan;  % Nan data from bad channels
    AnyArts(Interp,:)=0; Artifs.Epoch(Interp,:)=0; % Exlude those channels from trial rejection procedures 
    AnyBlinks = (Blinks.Epoch)>0;
    Exclude = AnyBlinks>0 | sum(AnyArts,1)>0;
    erp(:,:,Exclude) = nan; %erpr(:,:,Exclude)=nan;
    
    
    % Plot artifacts for each subject
    subplot(4,6,s); hold on;
    bar(1:128, sum(AnyArts(:,~AnyBlinks),2));
    numblinks = sum(AnyBlinks);
    numartsonly = sum(sum(AnyArts(:, ~AnyBlinks),1)>0);
    numtrialsdisp = sum(~Exclude);
    t=title([allsubj{s} ' ' num2str(numblinks) ' blinks ' num2str(numtrialsdisp) ' trials']); set(t,'fontsize', 14);
    
            %% Compile Single Trial artifact matrices

    if s==1,
        artifs.PreStim1=single(Artifs.PreStim1(:,1:size(erp,3))); 
        artifs.PreStim2=single(Artifs.PreStim2(:,1:size(erp,3)));
        artifs.Epoch=single(Artifs.Epoch(:,1:size(erp,3)));
        artifs.PostRTWind=single(Artifs.PostRTWind(:,1:size(erp,3))); 

        blinks.PreStim1=single(Blinks.PreStim1(:,1:size(erp,3)));
        blinks.PreStim2=single(Blinks.PreStim2(:,1:size(erp,3)));
        blinks.Epoch=single(Blinks.Epoch(:,1:size(erp,3)));
        blinks.PostRTWind=single(Blinks.PostRTWind(:,1:size(erp,3))); 
    else
        artifs.PreStim1=[artifs.PreStim1 single(Artifs.PreStim1(:,1:size(erp,3)))]; 
        artifs.PreStim2=[artifs.PreStim2 single(Artifs.PreStim2(:,1:size(erp,3)))];
        artifs.Epoch=[artifs.Epoch single(Artifs.Epoch(:,1:size(erp,3)))];
        artifs.PostRTWind=[artifs.PostRTWind single(Artifs.PostRTWind(:,1:size(erp,3)))]; 

        blinks.PreStim1=[blinks.PreStim1 single(Blinks.PreStim1(:,1:size(erp,3)))];
        blinks.PreStim2=[blinks.PreStim2 single(Blinks.PreStim2(:,1:size(erp,3)))];
        blinks.Epoch=[blinks.Epoch single(Blinks.Epoch(:,1:size(erp,3)))];
        blinks.PostRTWind=[blinks.PostRTWind single(Blinks.PostRTWind(:,1:size(erp,3)))]; 
    end
    
    clear Artifs Blinks

    %% Compile relevant characteristics for each single trial
    
    nn=nn+1; % A counter for all trials in the whole dataset. One matrix for all data from all subjects will be created.
    
    
    SubID(nn:nn+size(erp,3)-1)=single(s); Block(nn:nn+size(erp,3)-1) = single(blockid);
    Delay(nn:nn+size(erp,3)-1,:)=[single(delay); single(delay2)]'; Targ(nn:nn+size(erp,3)-1) = single(allTrig);
    Resp(nn:nn+size(erp,3)-1) = single(allrespLR); RT(nn:nn+size(erp,3)-1)=single(allRT);
    Acc(nn:nn+size(erp,3)-1)=single(Correct_Error); RTSamp(nn:nn+size(erp,3)-1)=single(allRTsamp);
    
    %% Compile FFT single trial matrices
    MotorChans = [33 34 48:56 61:64 106:117 122:124];
    FrontalChans = [75:78 84:91];
    OccChans = [4 8:10 14:24 27:31 37:39];
    
    fftwindlen = 205;
    Fres=fs/fftwindlen;
    slide = 25; % slide forward in time (ms)
    Freqs = 0:Fres:Fres*20;
    Steps(1) = round((length(slock_winds{1}) - fftwindlen)/slide);
    Steps(2) = round((length(slock_winds{2}) - fftwindlen)/slide);
    Steps(3) = round((length(slock_winds{3}) - fftwindlen)/slide);
    StepsR = round((length(rlock_winds) - fftwindlen)/slide);
    TimePoints = (round(slock_winds{3}(1)+((0:Steps(3))*slide)+(fftwindlen/2))) *(1000/512);
    TimePointsR = (round(rlock_winds(1)+((0:StepsR)*slide)+(fftwindlen/2))) *(1000/512);
    TimePointss = round(slock_winds{3}(1)+((0:Steps(3))*slide)+(fftwindlen/2));
    TimePointssR = round(rlock_winds(1)+((0:StepsR)*slide)+(fftwindlen/2));
    
    timefreq=[]; timefreqR=[];
    timefreq=nan(128,length(Freqs),length(TimePoints),size(erp,3), 'single');
    timefreqr=nan(128,length(Freqs),length(TimePointsR),size(erp,3), 'single');
    FFTMotor(1:length(MotorChans), 1:length(Freqs), 1:length(TimePointsR), nn:nn+size(erp,3)-1)=nan;
    FFTFrontal(1:length(FrontalChans), 1:length(Freqs), 1:length(TimePointsR), nn:nn+size(erp,3)-1)=nan;
    FFTOcc(1:length(OccChans), 1:length(Freqs), 1:length(TimePointsR), nn:nn+size(erp,3)-1)=nan;
    

    for n=1:size(erp,3)
            samp = find(t_dots{delay(n)}==delays(delay(n)));
            try motion_cohmo(:,nn+(n-1)) = Dots{n}(samp+ts_cohmo); DotMotion{nn+(n-1)}(:)=Dots{n}; catch; end

        if ~isnan(allRTsamp(n)) && allRTsamp(n)>(1000/(1000/512))
            if size(erp,2)>allRTsamp(n)+rlock_winds(end) % rlock window should not be out of bounds
                Rlockdata=erp(:,allRTsamp(n)+rlock_winds, n);
            else
                short=abs(size(erp,2)-(allRTsamp(n)+rlock_winds(end))); % In the event that RT+rlockwinds is out of bounds
                temprlock_winds=rlock_winds(1:(length(rlock_winds)-short));
                Rlockdata= erp(:,allRTsamp(n)+temprlock_winds, n);
            end
        end
        
        for step = 1:Steps(delay(n))+1
            thiswind = ismember(slock_winds{delay(n)},(slock_winds{delay(n)}+slide*(step-1)):(slock_winds{delay(n)}+slide*(step-1))+fftwindlen-1);
            fftdata=[]; fftdata = double(erp(:,thiswind, n));
            thisfft=[]; thisfft = single(abs(fft(fftdata, [],2))./(fftwindlen/2));
            FFTMotor(:,:,step, nn+(n-1))=thisfft(MotorChans,1:length(Freqs));
            FFTFrontal(:,:,step, nn+(n-1))=thisfft(FrontalChans,1:length(Freqs));
            FFTOcc(:,:,step, nn+(n-1))=thisfft(OccChans,1:length(Freqs));
            timefreq(:,:,step,n) = thisfft(:,1:length(Freqs));
            
            if ~isnan(allRT(n))
                StepsR=round(((size(Rlockdata,2)-205)/25)-.5);
                if step<StepsR+2
                    
                    thiswind=[]; thiswind = ismember(rlock_winds, (rlock_winds+slide*(step-1)):(rlock_winds+slide*(step-1))+fftwindlen-1);
                    fftdata=[]; fftdata = double(Rlockdata(:,thiswind));
                    thisfft=[]; thisfft = single(abs(fft(fftdata, [],2))./(fftwindlen/2));
                    FFTMotorR(:,:,step,nn+(n-1))=thisfft(MotorChans,1:length(Freqs));
                    FFTFrontalR(:,:,step, nn+(n-1))=thisfft(FrontalChans,1:length(Freqs));
                    FFTOccR(:,:,step, nn+(n-1))=thisfft(OccChans,1:length(Freqs));
                    timefreqR(:,:,step,n) = thisfft(:,1:length(Freqs));
                    
                end
            end
        end
    end
    
    timefreq=single(timefreq); timefreqR=single(timefreqR);
    nn=size(FFTMotor,4);
    
    clear erp;
    %%  ----------------------  Compile Matrices data averaged for each subject & condition ------------------- %%
    timefreq(:,:,:,Exclude) = nan;
        timefreqR(:,:,:,Exclude) = nan;
    
    for z=1:3
        for r=1:2
            numtrials(r,z,s) = sum(delay==z & allrespLR==RespCode(r) & ~Exclude);
            PlotData.TFreq(:,:,:,r,z,s)= single(squeeze(nanmean(timefreq(:, :,:, delay==z & allrespLR==RespCode(r)),4)));
            PlotData.TFreqr(:,:,:,r,z,s)= single(squeeze(nanmean(timefreqR(:, :,:, delay==z & allrespLR==RespCode(r)),4)));
            PlotData.RT(r,z,s) = nanmean(allRT(delay==z & allrespLR==RespCode(r)));
            PlotData.Acc(r,z,s) = sum(Acc(delay==z & allrespLR==RespCode(r))==1)/sum(delay==z & allrespLR==RespCode(r));
        end
    end
    
  
    numbins=3;
    for z=1:3
        for r=1:2
            rtlims(:,r,z,s) = prctile(allRT(delay==z & allrespLR==RespCode(r)), linspace(0, 100, numbins+1)) ;
            for bin=1:numbins
                PlotData.numtrialsbins(bin,r,z,s) = sum(allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s) & delay==z & allrespLR==RespCode(r) & allrespLR==RespCode(r));
                PlotData.RTBins(bin, r,z, s) = nanmean(allRT(allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s) & delay==z & allrespLR==RespCode(r) & allrespLR==RespCode(r))); % Get average RT within each RTbin
                PlotData.AccBins(bin,r,z,s) = sum(delay==z & allrespLR==RespCode(r) & allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s) & Correct_Error==1)/sum(delay==z & allrespLR==RespCode(r) & allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s));
                PlotData.TimeFreqbins(:,:,:,bin,r,z,s)= single(squeeze(nanmean(timefreq(:, :, :,delay==z & allrespLR==RespCode(r) & allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s)),4)));
                PlotData.TimeFreqrbins(:,:,:,bin,r,z,s)= single(squeeze(nanmean(timefreqR(:, :, :,delay==z & allrespLR==RespCode(r) & allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s)),4)));
            end
        end
    end
    
    %%
    clear erp allRT allrespLR delay delay2 sessid blockid Correct_Error allRTsamp
    FFTMotor=single(FFTMotor); FFTMotorR=single(FFTMotorR);
    FFTOcc=single(FFTOcc); FFTOccR=single(FFTOccR);
    FFTFrontal=single(FFTFrontal); FFTFrontalR=single(FFTFrontalR);
    disp(['subject ' allsubj{s} ' HPF ' Files{HPF} ' CSD ' num2str(spline)]) % display number of valid trials (ie RT greater than 75 but less than 1900
end

SubID=single(SubID); Block=single(Block); Delay=single(Delay); Targ=single(Targ);Resp=single(Resp); RT=single(RT); Acc=single(Acc); RTSamp=single(RTSamp);
%% Save data
if ismember(spline, 2:7)
    save([matfolder2 SubFold{HPF} 'SingleTrial_MotorFFT_' Files{HPF} '_CSD' num2str(spline)],'FFTMotor', 'FFTMotorR', 'MotorChans', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_FrontalFFT_' Files{HPF} '_CSD' num2str(spline)],'FFTFrontal', 'FFTFrontalR', 'FrontalChans','fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR', 'artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_OccFFT_' Files{HPF} '_CSD' num2str(spline)],'FFTOcc', 'FFTOccR', 'OccChans','fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR', 'artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([plotfolder, 'GrAv_TimeFreq' Files{HPF} '_CSD' num2str(spline),'.mat'], 'PlotData','numtrials', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','slock_wind', 'rlock_wind', 'delays', 'rtlims');
    save([matfolder2 SubFold{HPF} 'Dots'],'motion_cohmo', 'DotMotion', 'Matches','TrialMatch');
elseif spline==1
    save([matfolder2 SubFold{HPF} 'SingleTrial_MotorFFT_' Files{HPF} '_NoCSD'],'FFTMotor', 'FFTMotorR', 'MotorChans', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_FrontalFFT_' Files{HPF} '_NoCSD'],'FFTFrontal', 'FFTFrontalR', 'FrontalChans', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_OccFFT_' Files{HPF} '_NoCSD'],'FFTOcc', 'FFTOccR', 'OccChans','fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR', 'artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([plotfolder, 'GrAv_TimeFreq' Files{HPF} '_NoCSD.mat'], 'PlotData','numtrials', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','slock_wind', 'rlock_wind', 'delays', 'rtlims');
end


clear all
close all
