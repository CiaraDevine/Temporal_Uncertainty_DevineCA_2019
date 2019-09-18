clear all
close all

%% File and folder names

matfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\PreProcessedData\']; % Where to look for incoming data
SubFold = {'LPF_NoHPF\' 'LPF_HPF_0_05\' 'LPF_HPF_0_15\' 'LPF_HPF_0_25\' 'LPF_HPF_0_35\' 'LPF_HPF_0_45\'  'LPF_HPF_0_55\'  'LPF_HPF_0_65\'}; % Which subfolder to look in depending on which high pass filter being used

matfolder2 = ['C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\ERPMatrices\']; % Where to put outgoing data
Files = {'_NoHPF' '_HPF_0_05' '_HPF_0_15' '_HPF_0_25' '_HPF_0_35' '_HPF_0_45' '_HPF_0_55' '_HPF_0_65'}; % final part of file name changes depending on which HPF is being used.

plotfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\DotsPlots\'];
dotsfolder = 'C:\Users\devineca\Desktop\Dropbox\Ciara_DotsDelays\Data\DotsMats\';

%%



allsubj=  { 'GD'  'KM'  'OM'  'AF' 'SB' 'DH' 'SM' 'SD' 'HM' 'KW' 'DM' 'AE' 'ST' 'KL' 'SY' 'WR' 'JE' 'CM' 'EC' 'DPM' 'TS' 'SCM'}; % In order of when they were tested
blocks = {[1:5] [1:6] [1:5] [1:10] [1:10] [1:10] [1:10] [1:10] [1:2 4:10] [1:9] [1:9] [1:8] [1:8] [1:7] [1:8] [1:9] [1:8] [1:8] [1:8] [1:8] [1:8] [1:8]}; % For HM Block 3 got split into block 3 & 4 due to recording issue


fs = 512;
[B,A]=butter(4,6*2/fs); % creating filter for ERP waveforms
delays=[800 1200 1600];

%% Dots parameters
Fs = 100;

t_dots = {[10:10:2800],[10:10:3200],[10:10:3600]};
ts_dots={};
for delay = 1:3
    ts_dots{delay} = t_dots{delay}*Fs/1000;
end

rlockdotss = [-1*Fs:0*Fs];
rlockdots = rlockdotss*1000/Fs;
trs = [-0.750*Fs:0.125*Fs];
tr = trs*1000/Fs;



spline=4;
HPF=2;

f1=figure; hold on; st=suptitle('Artifacts'); set(st, 'fontsize',20, 'fontweight', 'bold');

%% Loop through subjects
SubID=[];Delay=[]; nn=0; Block=[]; RT=[]; Acc=[]; targ=[]; resp=[]; ERP=[]; ERPr=[];
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
    
    slock_winds=tts; slock_wind=tt;
    erp(erp==0)=NaN;
    
    
    %% Exclude artifact trials
    
    AnyBlinks = []; AnyArts=[]; Exclude=[]; ZERPEv=[]; ZERPRT=[];
    
    AnyArts = (Artifs.Epoch)>0;
    Interp=find(sum(AnyArts,2)>30); % Check for really bad single channels that should have been interpolated
    erp(Interp,:,:)=nan; erpr(Interp,:,:)=nan;  % Nan data from bad channels
    AnyArts(Interp,:)=0; Artifs.Epoch(Interp,:)=0; % Exlude those channels from trial rejection procedures 
    AnyBlinks = (Blinks.Epoch)>0;
    Exclude = AnyBlinks>0 | sum(AnyArts,1)>0;
    erp(:,:,Exclude) = nan; erpr(:,:,Exclude)=nan;
    
    
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
    
    %% Compile Single trial CPP matrices
    CPPch = [1 2 3 4 5 18 19 20 21 31 32 33 34 111 112 ];
    LRPch=[49:55 62:64 107:116 123 124];
    SPNch = [7:10 125:127 36:39 43:45];
    FCNch = [75:77 84:90];
    

    CPP(:,:,nn:nn+size(erp,3)-1)=single(erp(CPPch, 1:length(slock_wind{3}),:));% Channels x Time Points x Total Trials matrix. Retain only standard CPP electrodes (other wise computer might explode!).
    LRP(:,:,nn:nn+size(erp,3)-1)=single(erp(LRPch, 1:length(slock_wind{3}),:));
    SPN(:,:,nn:nn+size(erp,3)-1)=single(erp(SPNch, 1:length(slock_wind{3}),:));% Channels x Time Points x Total Trials matrix. Retain only standard CPP electrodes (other wise computer might explode!).
    FCN(:,:,nn:nn+size(erp,3)-1)=single(erp(FCNch, 1:length(slock_wind{3}),:));
    
   
    erpr=[]; DotMotionR(1:length(rlockdotss), nn:nn+size(erp,3)-1)=nan;
    for n=1:size(erp,3)
        rlocksamp=[];
        
        try 
            DotMotion{nn+(n-1)}(:)=Dots{n}; 
            rlocksamp = round(round(slock_wind{3}(allRTsamp(n)))/10);
            DotMotionR(:,nn+(n-1)) = Dots{n}(rlocksamp+rlockdotss); 
        catch
            DotMotion{nn+(n-1)}(:)=nan; 
        end
        
        if RT(nn+(n-1))>0
            erpr(:,:,n)=single(erp(:,allRTsamp(n)+rlock_winds, n));
            %RTSamp = find(slock_winds{3}==0) + delay2(n) + allRT(n);
            CPPr(:,:, nn+(n-1)) =  single(erpr(CPPch,:, n));
            LRPr(:,:, nn+(n-1)) =  single(erpr(LRPch,:, n));
            SPNr(:,:, nn+(n-1)) =  single(erpr(SPNch,:, n));
            FCNr(:,:, nn+(n-1)) =  single(erpr(FCNch,:, n));
        else
            erpr(:,:, n) = single(NaN(128,length(rlock_winds)));
            CPPr(:,:, nn+(n-1)) = single(NaN(length(CPPch),length(rlock_winds)));
            LRPr(:,:, nn+(n-1)) = single(NaN(length(LRPch),length(rlock_winds)));
            SPNr(:,:, nn+(n-1)) = single(NaN(length(SPNch),length(rlock_winds)));
            FCNr(:,:, nn+(n-1)) = single(NaN(length(FCNch),length(rlock_winds)));
        end
    end
    
    nn=size(CPP,3);
    
    
    
    %%  ----------------------  Compile Matrices data averaged for each subject & condition ------------------- %%


    %% Check for & Exclude outliers
    % Need to normalise EEG before computing z scores - baseline correct first 
    erprBL=erpr - repmat(nanmean(erp(:,slock_wind{3}>600 & slock_wind{3}<650,:),2), [1 size(erpr,2) 1]);
    ERPRT = squeeze(nanmean(erprBL(:, rlock_wind>-100 & rlock_wind<100, :),2));
    
    ZERPRT=nan(128,length(ERPRT));
    for z=1:3
        %ZERPEv(:,delay==z & ~Exclude) = squeeze(zscore(ERPEvOn(:,delay==z & ~Exclude),[],2));
        ZERPRT(:,delay==z & ~Exclude & allrespLR>0) = squeeze(zscore(ERPRT(:,delay==z & ~Exclude & allrespLR>0),[],2));
    end
    Outliers = ZERPEv>3 | ZERPEv<-3;
    OutlierIdx = repmat(permute(Outliers, [1 3 2]), [1, size(erp,2),1]);
    erp(OutlierIdx)=NaN; clear OutlierIdx
    OutlierIdx = repmat(permute(Outliers, [1 3 2]), [1, size(erpr,2),1]);
    erpr(OutlierIdx)=nan; clear OutlierIdx
    PlotData.Outliers(:,s) = sum(Outliers,2);
    
    
    %% Compile averaged ERP matrices
    for z=1:3
        PlotData.numtrials(z,s) = sum(delay==z & ~Exclude);
        PlotData.AvERP(:,:,z,s)= single(squeeze(nanmean(erp(:, :, delay==z),3)));
        PlotData.AvERPr(:,:,z,s)= single(squeeze(nanmean(erpr(:, :, delay==z),3)));
        PlotData.RT(z,s) = nanmean(allRT(delay==z & ~Exclude));
        PlotData.Acc(z,s) = sum(Correct_Error(delay==z)==1)/sum(delay==z);
    end
    

    clear rtlims;
    numbins=2;
    for z=1:3
        rtlims(:,z,s) = prctile(allRT(delay==z), linspace(0, 100, numbins+1)) ;
        for bin=1:numbins
            PlotData.RTBins1(bin, z, s) = nanmean(allRT(allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s) & delay==z)); % Get average RT within each RTbin
            PlotData.AccBins1(bin,z,s) = sum(delay(1,:)==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s) & Correct_Error==1)/sum(delay==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s));
            PlotData.CPPbins1(:,:,bin,z,s)= single(squeeze(nanmean(erp(:, :, delay==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s)),3)));
            PlotData.CPPrbins1(:,:,bin,z,s)= single(squeeze(nanmean(erpr(:, :, delay==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s)),3)));
            PlotData.NumtrialsRTBins1(bin,z,s)= sum(delay(1,:)==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s) & ~Exclude);
        end
    end
    PlotData.rtlims1(:,:,s)=rtlims(:,:,s);
    clear rtlims;
    
    numbins=3;
    for z=1:3
        rtlims(:,z,s) = prctile(allRT(delay==z), linspace(0, 100, numbins+1)) ;
        for bin=1:numbins
            PlotData.RTBins2(bin, z, s) = nanmean(allRT(allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s) & delay==z)); % Get average RT within each RTbin
            PlotData.AccBins2(bin,z,s) = sum(delay(1,:)==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s) & Correct_Error==1)/sum(delay==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s));
            PlotData.CPPbins2(:,:,bin,z,s)= single(squeeze(nanmean(erp(:, :, delay==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s)),3)));
            PlotData.CPPrbins2(:,:,bin,z,s)= single(squeeze(nanmean(erpr(:, :, delay==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s)),3)));
            PlotData.NumtrialsRTBins2(bin,z,s)= sum(delay(1,:)==z & allRT>rtlims(bin,z,s) & allRT<rtlims(bin+1,z,s) & ~Exclude);
        end
    end
    PlotData.rtlims2(:,:,s)=rtlims(:,:,s);
    
    PlotData.FastSlow500CPP(:,:,1,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT<250),3)));
    PlotData.FastSlow500CPP(:,:,2,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT>250),3)));
    PlotData.NumFastResp500(s) = sum(delay==3 & allRT<250 & ~Exclude);
    PlotData.FastSlow400CPP(:,:,1,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT<200),3)));
    PlotData.FastSlow400CPP(:,:,2,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT>200),3)));
    PlotData.NumFastResp400(s) = sum(delay==3 & allRT<200 & ~Exclude);
    PlotData.FastSlow350CPP(:,:,1,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT<175),3)));
    PlotData.FastSlow350CPP(:,:,2,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT>175),3)));
    PlotData.NumFastResp350(s) = sum(delay==3 & allRT<175 & ~Exclude);
    PlotData.FastSlow300CPP(:,:,1,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT<150),3)));
    PlotData.FastSlow300CPP(:,:,2,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT>150),3)));
    PlotData.NumFastResp300(s) = sum(delay==3 & allRT<150 & ~Exclude);
    PlotData.FastSlow250CPP(:,:,1,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT<125),3)));
    PlotData.FastSlow250CPP(:,:,2,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT>125),3)));
    PlotData.NumFastResp250(s) = sum(delay==3 & allRT<125 & ~Exclude);
    
    
    %%
    clear erp erpr allRT allrespLR delay delay2 sessid blockid Correct_Error allRTsamp
    disp(['subject ' allsubj{s} ' HPF ' Files{HPF} ' CSD ' num2str(spline)]) % display number of valid trials (ie RT greater than 75 but less than 1900
    disp([num2str(size(CPP,3)) ' ' num2str(size(DotMotion))]) % display number of valid trials (ie RT greater than 75 but less than 1900

end

CPP=single(CPP); CPPr=single(CPPr);
SubID=single(SubID); Block=single(Block); Delay=single(Delay); Targ=single(Targ);
Resp=single(Resp); RT=single(RT); Acc=single(Acc); RTSamp=single(RTSamp);

%% Save data
if ismember(spline, 2:7)
    save([matfolder2 SubFold{HPF} 'SingleTrial_CPP_' Files{HPF} '_CSD' num2str(spline)],'CPP','CPPr','CPPch','artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_LRP_' Files{HPF} '_CSD' num2str(spline)],'LRP','LRPr','LRPch','artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_SPN_' Files{HPF} '_CSD' num2str(spline)],'SPN','SPNr','SPNch','artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_FCN_' Files{HPF} '_CSD' num2str(spline)],'FCN','FCNr','FCNch','artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'Dots'], 'DotMotion', 'DotMotionR', 'Matches','TrialMatch', 'Matches');
    save([plotfolder, 'GrAv_ERP' Files{HPF} '_CSD' num2str(spline), '.mat'], 'PlotData', 'slock_wind', 'rlock_wind', 'delays', 'rtlims');
    clear PlotData;
elseif spline==1
    save([matfolder2 SubFold{HPF} 'SingleTrial_CPP_' Files{HPF} '_NoCSD_BL' num2str(BL)],'CPP','CPPr', 'CPPchans', 'artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_LRP_' Files{HPF} '_NoCSD_BL' num2str(BL)],'LRP','LRPr', 'LRPchans', 'artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_SPN_' Files{HPF} '_NoCSD_BL' num2str(BL)],'SPN','SPNr', 'SPNchans', 'artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([matfolder2 SubFold{HPF} 'SingleTrial_FCN_' Files{HPF} '_NoCSD_BL' num2str(BL)],'FCN','FCNr', 'LRPchans', 'artifs','blinks','SubID', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
    save([plotfolder, 'GrAv_ERP' Files{HPF} '_NoCSD.mat'], 'PlotData', 'slock_wind', 'rlock_wind', 'delays', 'rtlims');
    clear PlotData;
end

savefig([matfolder SubFold{HPF} 'Artifacts_CSD' num2str(spline) '_' num2str(Files{HPF})]);

clear all
close all
