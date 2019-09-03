clear all
close all

%% File and folder names
%path('C:\Users\devineca\Desktop\Learning_data\Analysis\Matlab_scripts\EEG\');
matfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\EpochData\']; % Where to look for incoming data
SubFold = {'LPF_NoHPF\' 'LPF_HPF_0_05\' 'LPF_HPF_0_15\' 'LPF_HPF_0_25\'}; % Which subfolder to look in depending on which high pass filter being used
matfolder2 = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\SingleTrialMatrices\']; % Where to put outgoing data
Files = {'_NoHPF' '_HPF_0_05' '_HPF_0_15' '_HPF_0_25'}; % final part of file name changes depending on which HPF is being used.

%%

allsubj =  { 'CJG'  'CR'  'AR' 'KD' 'TB'  'AOB'  'PM'  'CMG'  'SR'  'AB' 'GK' 'ED' 'SC' 'CB' 'CE' 'SH' 'RMC', 'JH', 'JR', 'ROC' 'KM' 'ROS'}; % In order of when they were tested %JR and KD excluded
sessions = {[1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:4] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5]  [1:5]};


fs = 512;
[B,A]=butter(4,6*2/fs); % creating filter for ERP waveforms
delays=[800 1200 1600];

RespCode=[12 13];
HPF=2;
spline=4; 


for k=[1 2]
    clear FFT*
    
    figure; hold on; st=suptitle('Artifacts'); set(st, 'fontsize',20, 'fontweight', 'bold');
    %% Loop through subjects
    FFTMotor=[]; FFTMotorR=[]; FFTOcc=[]; FFTOccR=[];FFTFrontal=[]; FFTFrontalR=[];
    SubID=[];Delay=[];  Block=[]; RT=[]; Acc=[]; Targ=[]; Resp=[]; Sessid=[]; RTSamp=[];
    nn=0;
    for s=1:length(allsubj)
        
        % Load in Data
        if spline ==1
            load([matfolder, SubFold{HPF}, allsubj{s}, 'GratingsDelaysERP_NoCSD', Files{HPF}]);
        else
            load([matfolder, SubFold{HPF}, allsubj{s}, 'GratingsDelaysERP_CSDSpl',num2str(spline), '_', Files{HPF}]);
        end
        allRTsamp=allRTsamp(1:size(erp,3));
        slock_winds = tts; slock_wind = tt; clear tt tts

        rlock_winds=-1*fs:.6*fs;
        rlock_wind=rlock_winds*(1000/fs);
        erp=erp(:,1:length(slock_winds{3}),:);
        erp(erp==0)=NaN;

        %% Exclude artifact trials
    
    AnyBlinks = []; AnyArts=[]; Exclude=[]; ZERPEv=[]; ZERPRT=[];
    
    AnyArts = Artifs.Epoch>0;
    for d=sessions{s}
    Interp=find(sum(AnyArts(:,sessid==d),2)>sum(sessid==d)/10); % Check for really bad single channels that should have been interpolated
    erp(Interp,:,sessid==d)=nan; %erpr(Interp,:,:)=nan;  % Nan data from bad channels
    AnyArts(Interp,sessid==d)=0; Artifs.Epoch(Interp,sessid==d)=0; % Exlude those channels from trial rejection procedures 
    end
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
        
        SubID(nn:nn+size(erp,3)-1)=single(s); 
        Block(nn:nn+size(erp,3)-1) = single(blockid); 
        Sessid(nn:nn+size(erp,3)-1)=single(sessid);
        Delay(nn:nn+size(erp,3)-1,:)=[single(delay); single(delay2)]'; 
        Targ(nn:nn+size(erp,3)-1) = single(allTrig);
        Resp(nn:nn+size(erp,3)-1) = single(allrespLR); 
        RT(nn:nn+size(erp,3)-1)=single(allRT);
        Acc(nn:nn+size(erp,3)-1)=single(Correct_Error); 
        RTSamp(nn:nn+size(erp,3)-1)=single(allRTsamp);
        
        
        Sessid=single(Sessid); SubID=single(SubID); Block=single(Block); Delay=single(Delay); 
        Targ=single(Targ); Resp=single(Resp); RT=single(RT); Acc=single(Acc); RTSamp=single(RTSamp);
        %% Loop through all trials
        
        %% Compile FFT single trial matrices
        MotorChans = [33 34 49:56 62:64 107:117 123 124];
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
        timefreqR=nan(128,length(Freqs),length(TimePointsR),size(erp,3), 'single');
        if k==1, FFTMotor(1:length(MotorChans), 1:length(Freqs), 1:length(TimePoints), nn:nn+size(erp,3)-1)=single(nan);
        elseif k==2, FFTFrontal(1:length(FrontalChans), 1:length(Freqs), 1:length(TimePoints), nn:nn+size(erp,3)-1)=single(nan);
        elseif k==3, FFTOcc(1:length(OccChans), 1:length(Freqs), 1:length(TimePoints), nn:nn+size(erp,3)-1)=single(nan);
        end
        
        for n=1:size(erp,3)
            if ~isnan(allRTsamp(n)) && allRTsamp(n)>(1000/(1000/512))
                if size(erp,2)>allRTsamp(n)+rlock_winds(end) % rlock window should not be out of bounds 
                Rlockdata=erp(:,allRTsamp(n)+rlock_winds, n);      
                else
                short=abs(size(erp,2)-(allRTsamp(n)+rlock_winds(end))); % In the event that RT+rlockwinds is out of bounds 
                temprlock_winds=rlock_winds(1:(length(rlock_winds)-short));
                Rlockdata= erp(:,allRTsamp(n)+temprlock_winds, n);
                end
            else 
                fluff=fluff;
            end
                        
            for step = 1:Steps(delay(n))+1
                thiswind = ismember(slock_winds{delay(n)},(slock_winds{delay(n)}+slide*(step-1)):(slock_winds{delay(n)}+slide*(step-1))+fftwindlen-1);
                fftdata=[]; fftdata = double(erp(:,thiswind, n));
                thisfft=[]; thisfft = single(abs(fft(fftdata, [],2))./(fftwindlen/2));
                if k==1, FFTMotor(:,:,step, nn+(n-1))=thisfft(MotorChans,1:length(Freqs));
                elseif k==2, FFTFrontal(:,:,step, nn+(n-1))=thisfft(FrontalChans,1:length(Freqs));
                elseif k==3, FFTOcc(:,:,step, nn+(n-1))=thisfft(OccChans,1:length(Freqs));
                end
                timefreq(:,:,step,n) = thisfft(:,1:length(Freqs));
                
                if ~isnan(allRT(n))
                    StepsR=round(((size(Rlockdata,2)-205)/25)-.5);
                    if step<StepsR+2
                        
                        thiswind=[]; thiswind = ismember(rlock_winds, (rlock_winds+slide*(step-1)):(rlock_winds+slide*(step-1))+fftwindlen-1);
                        fftdata=[]; fftdata = double(Rlockdata(:,thiswind));
                        thisfft=[]; thisfft = single(abs(fft(fftdata, [],2))./(fftwindlen/2));
                        if k==1, FFTMotorR(:,:,step,nn+(n-1))=thisfft(MotorChans,1:length(Freqs));
                        elseif k==2,  FFTFrontalR(:,:,step, nn+(n-1))=thisfft(FrontalChans,1:length(Freqs));
                        elseif k==3,  FFTOccR(:,:,step, nn+(n-1))=thisfft(OccChans,1:length(Freqs));
                        end
                        timefreqR(:,:,step,n) = thisfft(:,1:length(Freqs));
                        
                    end
                end
            end
        end
        
        if k==1, nn=size(FFTMotor,4);
        elseif k==2, nn=size(FFTFrontal,4);
        elseif k==3, nn=size(FFTOcc,4);
        end
        
        clear erp;
        
        %%  ----------------------  Compile Matrices data averaged for each subject & condition ------------------- %%
        
        
        timefreq(:,:,:,Exclude) = nan;
        timefreqR(:,:,:,Exclude) = nan;
        
        % Compile averaged ERP matrices
        delay_Npl=[]; delay_Npl=[delay(2:length(delay)) nan]; block_Npl=[]; block_Npl=[blockid(2:length(blockid)) nan]; 

        for z=1:3
            for r=1:2
                PlotData.numtrials(r,z,s) = sum(delay==z & allrespLR==RespCode(r) & ~Exclude);
                PlotData.TFreq(:,:,:,r,z,s)= single(squeeze(nanmean(timefreq(:, :,:, delay==z & allrespLR==RespCode(r)),4)));
                PlotData.TFreqr(:,:,:,r,z,s)= single(squeeze(nanmean(timefreqR(:, :,:, delay==z & allrespLR==RespCode(r)),4)));
                PlotData.RT(r,z,s) = nanmean(RT(delay==z & allrespLR==RespCode(r)));
                PlotData.Acc(r,z,s) = sum(Correct_Error(delay==z & allrespLR==RespCode(r))==1)/sum(delay==z & allrespLR==RespCode(r));
                
                PlotData2.numtrials(r,z,s) = sum(delay==z & allTrig==r+7 & ~Exclude);
                PlotData2.TFreq(:,:,:,r,z,s)= single(squeeze(nanmean(timefreq(:, :,:, delay==z & allTrig==r+7),4)));
                PlotData2.TFreqr(:,:,:,r,z,s)= single(squeeze(nanmean(timefreqR(:, :,:, delay==z & allTrig==r+7),4)));
                PlotData2.RT(r,z,s) = nanmean(RT(delay==z & allTrig==r+7));
                PlotData2.Acc(r,z,s) = sum(Correct_Error(delay==z & allTrig==r+7)==1)/sum(delay==z & allTrig==r+7);
                for zz=1:3
                    PlotData3.numtrials(r,zz,z,s) = sum(delay==z & delay_Npl==zz & blockid==block_Npl &  allrespLR==RespCode(r) & ~Exclude);
                    PlotData3.TFreq(:,:,:,r,zz,z,s)= single(squeeze(nanmean(timefreq(:, :,:, delay==z & delay_Npl==zz & blockid==block_Npl & allrespLR==RespCode(r)),4)));
                    PlotData3.TFreqr(:,:,:,r,zz,z,s)= single(squeeze(nanmean(timefreqR(:, :,:, delay==z & delay_Npl==zz & blockid==block_Npl & allrespLR==RespCode(r)),4)));
                    PlotData3.RT(r,z,s) = nanmean(RT(delay==z & delay_Npl==zz & blockid==block_Npl & allrespLR==RespCode(r)));
                    PlotData3.Acc(r,z,s) = sum(Correct_Error(delay==z & delay_Npl==zz & blockid==block_Npl & allrespLR==RespCode(r))==1)/sum(delay==z & delay_Npl==zz & blockid==block_Npl & allrespLR==RespCode(r));
                end
                
                for d=1:5
                    PlotTraining1.numtrials(r,z,d,s) = sum(sessid==d & delay==z & allrespLR==RespCode(r) & ~Exclude);
                    PlotTraining1.TFreq(:,:,:,r,z,d,s)= single(squeeze(nanmean(timefreq(:, :,:, sessid==d & delay==z & allrespLR==RespCode(r)),4)));
                    PlotTraining1.TFreqr(:,:,:,r,z,d,s)= single(squeeze(nanmean(timefreqR(:, :,:, sessid==d & delay==z & allrespLR==RespCode(r)),4)));
                    PlotTraining1.RT(r,z,d,s) = nanmean(RT(sessid==d & delay==z & allrespLR==RespCode(r)));
                    PlotTraining1.Acc(r,z,d,s) = sum(Correct_Error(sessid==d & delay==z & allrespLR==RespCode(r))==1)/sum(sessid==d & delay==z & allrespLR==RespCode(r));
                    
                    PlotTraining2.numtrials(r,z,d,s) = sum(sessid==d & delay==z & allTrig==r+7 & ~Exclude);
                    PlotTraining2.TFreq(:,:,:,r,z,d,s)= single(squeeze(nanmean(timefreq(:, :,:, sessid==d & delay==z & allTrig==r+7),4)));
                    PlotTraining2.TFreqr(:,:,:,r,z,d,s)= single(squeeze(nanmean(timefreqR(:, :,:, sessid==d & delay==z & allTrig==r+7),4)));
                    PlotTraining2.RT(r,z,d,s) = nanmean(RT(sessid==d & delay==z & allTrig==r+7));
                    PlotTraining2.Acc(r,z,d,s) = sum(Correct_Error(sessid==d & delay==z & allTrig==r+7)==1)/sum(sessid==d & delay==z & allTrig==r+7);
                    
                end
            end
        end
        
        PlotData.TFreqFast(:,:,:,1,1,s) = single(squeeze(nanmean(timefreq(:, :,:, delay==z & allrespLR==RespCode(1) & allRT<250),4)));
        PlotData.TFreqFast(:,:,:,2,1,s) = single(squeeze(nanmean(timefreq(:, :,:, delay==z & allrespLR==RespCode(2) & allRT<250),4)));
        PlotData.TFreqFast(:,:,:,1,2,s) = single(squeeze(nanmean(timefreq(:, :,:, delay==z & allrespLR==RespCode(1) & allRT>250),4)));
        PlotData.TFreqFast(:,:,:,2,2,s) = single(squeeze(nanmean(timefreq(:, :,:, delay==z & allrespLR==RespCode(2) & allRT>250),4)));
        PlotData2.TFreqFast(:,:,:,1,1,s) = single(squeeze(nanmean(timefreq(:, :,:, delay==z & allTrig==8 & allRT<250),4)));
        PlotData2.TFreqFast(:,:,:,2,1,s) = single(squeeze(nanmean(timefreq(:, :,:, delay==z & allTrig==9 & allRT<250),4)));
        PlotData2.TFreqFast(:,:,:,1,2,s) = single(squeeze(nanmean(timefreq(:, :,:, delay==z & allTrig==8 & allRT>250),4)));
        PlotData2.TFreqFast(:,:,:,2,2,s) = single(squeeze(nanmean(timefreq(:, :,:, delay==z & allTrig==9 & allRT>250),4)));
        
        numbins=6;
        for z=1:3
            for r=1:2
                rtlims(:,r,z,s) = prctile(allRT(delay==z & allrespLR==RespCode(r)), linspace(0, 100, numbins+1)) ;
                rtlims2(:,r,z,s) = prctile(allRT(delay==z & allTrig==r+7), linspace(0, 100, numbins+1)) ;
                for bin=1:numbins
                    PlotData.numtrialsbins(bin,r,z,s) = sum(allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s) & delay==z & allrespLR==RespCode(r) & allrespLR==RespCode(r));
                    PlotData.RTBins(bin, r,z, s) = nanmean(allRT(allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s) & delay==z & allrespLR==RespCode(r) & allrespLR==RespCode(r))); % Get average RT within each RTbin
                    PlotData.AccBins(bin,r,z,s) = sum(delay==z & allrespLR==RespCode(r) & allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s) & Correct_Error==1)/sum(delay==z & allrespLR==RespCode(r) & allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s));
                    PlotData.TimeFreqbins(:,:,:,bin,r,z,s)= single(squeeze(nanmean(timefreq(:, :, :,delay==z & allrespLR==RespCode(r) & allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s)),4)));
                    PlotData.TimeFreqrbins(:,:,:,bin,r,z,s)= single(squeeze(nanmean(timefreqR(:, :, :,delay==z & allrespLR==RespCode(r) & allRT>rtlims(bin,r,z,s) & allRT<rtlims(bin+1,r,z,s)),4)));
                    
                    PlotData2.numtrialsbins(bin,r,z,s) = sum(allRT>rtlims2(bin,r,z,s) & allRT<rtlims2(bin+1,r,z,s) & delay==z & allTrig==r+7 & allrespLR==RespCode(r));
                    PlotData2.RTBins(bin, r,z, s) = nanmean(allRT(allRT>rtlims2(bin,r,z,s) & allRT<rtlims2(bin+1,r,z,s) & delay==z & allTrig==r+7)); % Get average RT within each RTbin
                    PlotData2.AccBins(bin,r,z,s) = sum(delay==z & allTrig==r+7 & allRT>rtlims2(bin,r,z,s) & allRT<rtlims2(bin+1,r,z,s) & Correct_Error==1)/sum(delay==z & allTrig==r+7 & allRT>rtlims2(bin,r,z,s) & allRT<rtlims2(bin+1,r,z,s));
                    PlotData2.TimeFreqbins(:,:,:,bin,r,z,s)= single(squeeze(nanmean(timefreq(:, :, :,delay==z & allTrig==r+7 & allRT>rtlims2(bin,r,z,s) & allRT<rtlims2(bin+1,r,z,s)),4)));
                    PlotData2.TimeFreqrbins(:,:,:,bin,r,z,s)= single(squeeze(nanmean(timefreqR(:, :, :,delay==z & allTrig==r+7 & allRT>rtlims2(bin,r,z,s) & allRT<rtlims2(bin+1,r,z,s)),4)));
                end
            end
        end
        PlotData.rtlims=rtlims;
        PlotData2.rtlims=rtlims2;
        
        
        %%
        
        
        if k==1, FFTMotor=single(FFTMotor); FFTMotorR=single(FFTMotorR);
        elseif k==2, FFTFrontal=single(FFTFrontal); FFTFrontalR=single(FFTFrontalR);
        elseif k==3, FFTOcc=single(FFTOcc); FFTOccR=single(FFTOccR);
        end
        
        clear erp* erpr* ERPRT ZERP* Exclude AnyArts AnyBlinks allRT allTrig allRTsamp allrespLR delay delay2 sessid blockid Correct_Error allRTsamp FixLen
        disp(['subject ' allsubj{s} ' HPF ' Files{HPF} ' CSD ' num2str(spline)]) % display number of valid trials (ie RT greater than 75 but less than 1900
        
        
    end

    Sessid=single(Sessid); SubID=single(SubID); Block=single(Block); Delay=single(Delay); Targ=single(Targ);
    Resp=single(Resp); RT=single(RT); Acc=single(Acc); RTSamp=single(RTSamp);

    %% Save data
    if ismember(spline, 2:7)
        if k==1, save([matfolder2 SubFold{HPF} 'SingleTrial_MotorFFT_' Files{HPF} '_CSD' num2str(spline)],'FFTMotor', 'FFTMotorR', 'MotorChans', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','artifs','blinks','SubID', 'Block','Sessid','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
        elseif k==2, save([matfolder2 SubFold{HPF} 'SingleTrial_FrontalFFT_' Files{HPF} '_CSD' num2str(spline)],'FFTFrontal', 'FFTFrontalR', 'FrontalChans','fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR', 'artifs','blinks','SubID', 'Block','Sessid','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
        elseif k==3, save([matfolder2 SubFold{HPF} 'SingleTrial_OccFFT_' Files{HPF} '_CSD' num2str(spline)],'FFTOcc', 'FFTOccR', 'OccChans','fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR', 'artifs','blinks','SubID', 'Block','Sessid','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
        end
        save([matfolder2, 'GrAv_TimeFreq' Files{HPF} '_CSD' num2str(spline) '.mat'], 'PlotData','PlotData2','PlotData3','PlotTraining1','PlotTraining2', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','slock_wind', 'rlock_wind', 'delays', 'rtlims');
    elseif spline==1
        if k==1, save([matfolder2 SubFold{HPF} 'SingleTrial_MotorFFT_' Files{HPF} '_NoCSD'],'FFTMotor', 'FFTMotorR', 'MotorChans', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','artifs','blinks','SubID', 'Block','Sessid','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
        elseif k==2, save([matfolder2 SubFold{HPF} 'SingleTrial_FrontalFFT_' Files{HPF} '_NoCSD'],'FFTFrontal', 'FFTFrontalR', 'FrontalChans', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','artifs','blinks','SubID', 'Block','Sessid','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
        elseif k==3, save([matfolder2 SubFold{HPF} 'SingleTrial_OccFFT_' Files{HPF} '_NoCSD'],'FFTOcc', 'FFTOccR', 'OccChans','fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR', 'artifs','blinks','SubID', 'Block','Sessid','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
        end
        save([matfolder2, 'GrAv_TimeFreq' Files{HPF} '_NoCSD.mat'], 'PlotData','PlotData2','PlotData3','PlotTraining1','PlotTraining2', 'fftwindlen', 'Fres', 'Freqs','TimePoint*','Steps', 'StepsR','slock_wind', 'rlock_wind', 'delays', 'rtlims');
    end
    
    savefig([matfolder2 SubFold{HPF} 'Artifacts_CSD' num2str(spline) '_' num2str(Files{HPF})]);
    

end