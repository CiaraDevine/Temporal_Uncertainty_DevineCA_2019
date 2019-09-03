clear all
close all

%% File and folder names
%path('C:\Users\devineca\Desktop\Learning_data\Analysis\Matlab_scripts\EEG\');
matfolder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\EpochedData\']; % Where to look for incoming data
SubFold = {'LPF_NoHPF\' 'LPF_HPF_0_05\' 'LPF_HPF_0_15\' 'LPF_HPF_0_25\'}; % Which subfolder to look in depending on which high pass filter being used
matfolder2 = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\SingleTrialMatrices\']; % Where to put outgoing data
Files = {'_NoHPF' '_HPF_0_05' '_HPF_0_15' '_HPF_0_25'}; % final part of file name changes depending on which HPF is being used.

%%

allsubj =  { 'CJG'  'CR'  'AR' 'KD' 'TB'  'AOB'  'PM'  'CMG'  'SR'  'AB' 'GK' 'ED' 'SC' 'CB' 'CE' 'SH' 'RMC', 'JH', 'JR', 'ROC' 'KM' 'ROS'}; % In order of when they were tested %JR and KD excluded
sessions = {[1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:4] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5]  [1:5]};


fs = 512;
[B,A]=butter(4,6*2/fs); % creating filter for ERP waveforms
delays=[800 1200 1600];

spline=4;
HPF=2;
RTwin=[-50 50];
BL=[550 650];

for k=2%[4 2 3 1]
    
    figure; hold on; st=suptitle('Artifacts'); set(st, 'fontsize',20, 'fontweight', 'bold');
    %% Loop through subjects
    SubID=[];Delay=[]; nn=0; Block=[]; RT=[]; Acc=[]; Targ=[]; Resp=[]; CPP=[]; CPPr=[]; Sessid=[]; RTSamp=[]; LRP=[]; LRPr=[]; SPN=[]; SPNr=[]; FCN=[];  FCNr=[];
    for s=1:length(allsubj)
        
        % Load in Data
        
        load([matfolder, SubFold{HPF}, allsubj{s}, 'GratingsDelaysERP_CSDSpl',num2str(spline), '_', Files{HPF}]);
        
        slock_winds = tts; slock_wind = tt; clear tt tts
        rlock_winds = round(-1*fs: .6*fs);rlock_wind = rlock_winds*(1000/512);
        erp=erp(:,1:length(slock_winds{3}),:);
        
        erp(erp==0)=NaN;
        %%
        allRTsamp = allRTsamp(1:length(blockid)); Correct_Error=single(Correct_Error);
        blockid=single(blockid); sessid=single(sessid); allRT=single(allRT);
        delay=single(delay); delay2=single(delay2); allRTsamp=single(allRTsamp);
        allrespLR=single(allrespLR); allTrig=single(allTrig);
        
        
        

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
        %artifs.PreStim1=single(Artifs.PreStim1(:,1:size(erp,3))); %artifs.PreStim2=single(Artifs.PreStim2(:,1:size(erp,3)));
        artifs.Epoch=single(Artifs.Epoch(:,1:size(erp,3))); %artifs.PostRTWind=single(Artifs.PostRTWind(:,1:size(erp,3))); 

        %blinks.PreStim1=single(Blinks.PreStim1(:,1:size(erp,3))); %blinks.PreStim2=single(Blinks.PreStim2(:,1:size(erp,3)));
        blinks.Epoch=single(Blinks.Epoch(:,1:size(erp,3))); %blinks.PostRTWind=single(Blinks.PostRTWind(:,1:size(erp,3))); 
    else
%         artifs.PreStim1=[artifs.PreStim1 single(Artifs.PreStim1(:,1:size(erp,3)))]; artifs.PreStim2=[artifs.PreStim2 single(Artifs.PreStim2(:,1:size(erp,3)))];
        artifs.Epoch=[artifs.Epoch single(Artifs.Epoch(:,1:size(erp,3)))]; %artifs.PostRTWind=[artifs.PostRTWind single(Artifs.PostRTWind(:,1:size(erp,3)))]; 

        %blinks.PreStim1=[blinks.PreStim1 single(Blinks.PreStim1(:,1:size(erp,3)))]; %blinks.PreStim2=[blinks.PreStim2 single(Blinks.PreStim2(:,1:size(erp,3)))];
        blinks.Epoch=[blinks.Epoch single(Blinks.Epoch(:,1:size(erp,3)))];% blinks.PostRTWind=[blinks.PostRTWind single(Blinks.PostRTWind(:,1:size(erp,3)))]; 
    end
    
clear Artifs Blinks
    
        
        %% Compile relevant characteristics for each single trial
        
        nn=nn+1; % A counter for all trials in the whole dataset. One matrix for all data from all subjects will be created.
        
        SubID(nn:nn+size(erp,3)-1)=single(s); Block(nn:nn+size(erp,3)-1) = single(blockid); Sessid(nn:nn+size(erp,3)-1)=single(sessid);
        Delay(nn:nn+size(erp,3)-1,:)=[single(delay); single(delay2)]'; Targ(nn:nn+size(erp,3)-1) = single(allTrig);
        Resp(nn:nn+size(erp,3)-1) = single(allrespLR); RT(nn:nn+size(erp,3)-1)=single(allRT)*(1000/512);
        Acc(nn:nn+size(erp,3)-1)=single(Correct_Error); RTSamp(nn:nn+size(erp,3)-1)=single(allRTsamp);

        % Reduce the memory taken by each matrix
        Sessid=single(Sessid); SubID=single(SubID); Block=single(Block); Delay=single(Delay); Targ=single(Targ);
        Resp=single(Resp); RT=single(RT); Acc=single(Acc); RTSamp=single(RTSamp);

        delayz=[800 1200 1600];
        %% Compile Single trial CPP matrices
        CPPch = [1:5 18 19 20 21 31 32 34 112];
        SPNch = [113 124:127 6:11 34:39 43:45 51];
        LRPch=[49:55 62:64 107:116 123 124];
        FCNch = [75:78 82:87 88:91];
        
        if k==1
            CPP(:,:,nn:nn+size(erp,3)-1)=single(erp(CPPch, 1:length(slock_wind{3}),:));% Channels x Time Points x Total Trials matrix. Retain only standard CPP electrodes (other wise computer might explode!).
        elseif k==2 
               erptemp=[]; 
               erptemp=nan(128, length(-51:409), size(erp,3));
            for z=1:3, 
                win=slock_wind{z}>delayz(z)-100 & slock_wind{z}<delays(z)+800; 
                erptemp(:,1:sum(win),delay==z)=erp(:, win, delay==z);
            end
            SPN(:,:,nn:nn+size(erp,3)-1)=single(erptemp(SPNch, :,:));% Channels x Time Points x Total Trials matrix. Retain only standard CPP electrodes (other wise computer might explode!).
        elseif k==3
            LRP(:,:,nn:nn+size(erp,3)-1)=single(erp(LRPch, 1:length(slock_wind{3}),:));
        elseif k==4 
            FCN(:,:,nn:nn+size(erp,3)-1)=single(erp(FCNch, 1:length(slock_wind{3}),:));
        end
        
        if k==1
            CPP=single(CPP);
        elseif k==2
            SPN=single(SPN);
        elseif k==3
            LRP=single(LRP);
        elseif k==4
            FCN=single(FCN);
        end
        
        
        erpr=nan(size(erp,1),length(rlock_wind),size(erp,3));
        for n=1:size(erp,3)
            if RT(nn+(n-1))>0
                if size(erp,2)>allRTsamp(n)+rlock_winds(end)
                    erpr(:,:,n)=single(erp(:,allRTsamp(n)+rlock_winds, n));
                else 
                   short=abs(size(erp,2)-(allRTsamp(n)+rlock_winds(end))); 
                   temprlock_winds=rlock_winds(1:(length(rlock_winds)-short));
                   erpr(:,1:(size(erpr,2)-short),n)=single(erp(:,allRTsamp(n)+temprlock_winds, n)); 
                end
                
                %RTSamp = find(slock_winds{3}==0) + delay2(n) + allRT(n);
                if k==1
                    CPPr(:,:, nn+(n-1)) =  single(erpr(CPPch,:, n));
                elseif k==2 
                    SPNr(:,:, nn+(n-1)) =  single(erpr(SPNch,:, n));
                elseif k==3
                    LRPr(:,:, nn+(n-1)) =  single(erpr(LRPch,:, n));
                elseif k==4
                    FCNr(:,:, nn+(n-1)) =  single(erpr(FCNch,:, n));
                end
            else
                erpr(:,:, n) = single(NaN(128,length(rlock_winds)));
                if k==1
                    CPPr(:,:, nn+(n-1)) = single(NaN(length(CPPch),length(rlock_winds)));
                elseif k==2    
                    SPNr(:,:, nn+(n-1)) = single(NaN(length(SPNch),length(rlock_winds)));
                elseif k==3
                    LRPr(:,:, nn+(n-1)) = single(NaN(length(LRPch),length(rlock_winds)));
                elseif k==4   
                    FCNr(:,:, nn+(n-1)) = single(NaN(length(FCNch),length(rlock_winds)));
                end
            end
        end
        erpr=single(erpr); erp=single(erp);
        if k==1
            CPPr=single(CPPr); nn=size(CPP,3);
        elseif k==2
            SPNr=single(SPNr);  nn=size(SPN,3);
        elseif k==3
            LRPr=single(LRPr); nn=size(LRP,3);
        elseif k==4
            FCNr=single(FCNr); nn=size(FCN,3);
        end
        
        
        
        
        
        %%  ----------------------  Compile Matrices data averaged for each subject & condition ------------------- %%
        

        %% Check for Outliers
        erprBL=[]; erprBL=erpr - repmat(nanmean(erp(:,slock_wind{3}>BL(1) & slock_wind{3}<BL(2),:),2), [1 size(erpr,2) 1]);
        ERPRT=[];ERPRT = squeeze(nanmean(erprBL(:, rlock_wind>RTwin(1) & rlock_wind<RTwin(2), :),2));
        
        ZERPRT=[]; ZERPRT=nan(128,length(ERPRT));
        checknans=[]; checknans=any(isnan(ERPRT));
        for d=sessions{s}
        for z=1:3
            idx=[]; idx=sessid==d & delay==z & ~Exclude & ~checknans;
            ZERPRT(:,idx) = squeeze(zscore(ERPRT(:,idx),[],2));
            Outlierstemp=[]; Outlierstemp= ZERPRT(:,idx)<-3 | ZERPRT(:,idx)>3;
            PlotData.Outliers(:,d,z,s) = sum(Outlierstemp,2);
        end
        end
        Outliers=[]; Outliers = ZERPRT<-3 | ZERPRT>3;
        OutlierIdx = repmat(permute(Outliers, [1 3 2]), [1, size(erp,2),1]);
        erp(OutlierIdx)=NaN; clear OutlierIdx
        OutlierIdx = repmat(permute(Outliers, [1 3 2]), [1, size(erpr,2),1]);
        erpr(OutlierIdx)=nan; clear OutlierIdx
         
        %% Compile averaged ERP matrices
        
        delay_Nmin=[]; block_Nmin=[]; 
        delay_Nmin=[nan delay(1:length(delay)-1)]; block_Nmin=[nan blockid(1:length(blockid)-1)]; 
        for z=1:3
            for zz=1:3
            PlotData.numtrials(zz,z,s) = sum(delay==z & delay_Nmin==zz & ~Exclude);
            PlotData.ERP(:,:,zz,z,s)= single(squeeze(nanmean(erp(:, :, delay==z & delay_Nmin==zz),3)));
            PlotData.ERPr(:,:,zz,z,s)= single(squeeze(nanmean(erpr(:, :, delay==z & delay_Nmin==zz),3)));
            PlotData.RT(zz,z,s) = nanmean(allRT(delay==z & delay_Nmin==zz & ~Exclude));
            PlotData.Acc(zz,z,s) = sum(Correct_Error(delay==z & delay_Nmin==zz)==1)/sum(delay==z & delay_Nmin==zz);
            for d=1:5
                PlotDataTraining.numtrials(d,zz,z,s) = sum(sessid==d & delay==z & delay_Nmin==zz & ~Exclude);
                PlotDataTraining.ERP(:,:,d,zz,z,s)= single(squeeze(nanmean(erp(:, :, sessid==d & delay==z & delay_Nmin==zz),3)));
                PlotDataTraining.ERPr(:,:,d,zz,z,s)= single(squeeze(nanmean(erpr(:, :, sessid==d & delay==z & delay_Nmin==zz),3)));
                PlotDataTraining.RT(d,zz,z,s) = nanmean(allRT(sessid==d & delay==z & delay_Nmin==zz));
                PlotDataTraining.Acc(d,zz,z,s) = sum(Correct_Error(sessid==d & delay==z & delay_Nmin==zz)==1)/sum(sessid==d & delay==z & delay_Nmin==zz);
            end
            end
        end
        
        
        for f=1; 
%         PlotData.FastSlowCPP500(:,:,1,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT<250),3)));
%         PlotData.FastSlowCPP500(:,:,2,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT>250),3)));
%         PlotData.NumFastCPP500(s) = sum(delay==3 & allRT<250 & ~Exclude);
% 
%         PlotData.FastSlowCPP250(:,:,1,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT<125),3)));
%         PlotData.FastSlowCPP250(:,:,2,s) = single(squeeze(nanmean(erp(:, :, delay==3 & allRT>125),3)));
%         PlotData.NumFastCPP250(s) = sum(delay==3 & allRT<125 & ~Exclude);
%         for d=1:5
%             PlotDataTraining.FastCPP(:,:,d,1,s) = single(squeeze(nanmean(erp(:, :, delay==3 & sessid==d & allRT<250),3)));
%             PlotDataTraining.FastCPP(:,:,d,2,s) = single(squeeze(nanmean(erp(:, :, delay==3 & sessid==d & allRT>250),3)));
%             PlotDataTraining.NumFastCPP(s) = sum(delay==3  & sessid==d & allRT<250 & ~Exclude);
%         end
%         
%         numbins1=6;
%         for z=1:3
%             rtlims1(:,z,s) = prctile(allRT(delay==z), linspace(0, 100, numbins1+1)) ;
%             for bin=1:numbins1
%                 PlotData.numtrialsBins(bin,z,s) = sum(allRT>rtlims1(bin,z,s) & allRT<rtlims1(bin+1,z,s) & delay==z & ~Exclude);
%                 PlotData.RTBins(bin, z, s) = nanmean(allRT(allRT>rtlims1(bin,z,s) & allRT<rtlims1(bin+1,z,s) & delay==z)); % Get average RT within each RTbin
%                 PlotData.AccBins(bin,z,s) = sum(delay(1,:)==z & allRT>rtlims1(bin,z,s) & allRT<rtlims1(bin+1,z,s) & Correct_Error==1)/sum(delay==z & allRT>rtlims1(bin,z,s) & allRT<rtlims1(bin+1,z,s));
%                 PlotData.ERPbins(:,:,bin,z,s)= single(squeeze(nanmean(erp(:, :, delay==z & allRT>rtlims1(bin,z,s) & allRT<rtlims1(bin+1,z,s)),3)));
%                 PlotData.ERPrbins(:,:,bin,z,s)= single(squeeze(nanmean(erpr(:, :, delay==z & allRT>rtlims1(bin,z,s) & allRT<rtlims1(bin+1,z,s)),3)));
%             end
%         end
%  
%         
%         numbins3=10;
%         for z=1:3
%             rtlims3(:,z,s) = prctile(allRT(delay==z), linspace(0, 100, numbins3+1)) ;
%             for bin=1:numbins3
%                 PlotData.numtrialsBins2(bin,z,s) = sum(allRT>rtlims3(bin,z,s) & allRT<rtlims3(bin+1,z,s) & delay==z & ~Exclude);
%                 PlotData.RTBins2(bin, z, s) = nanmean(allRT(allRT>rtlims3(bin,z,s) & allRT<rtlims3(bin+1,z,s) & delay==z)); % Get average RT within each RTbin
%                 PlotData.AccBins2(bin,z,s) = sum(delay(1,:)==z & allRT>rtlims3(bin,z,s) & allRT<rtlims3(bin+1,z,s) & Correct_Error==1)/sum(delay==z & allRT>rtlims3(bin,z,s) & allRT<rtlims3(bin+1,z,s));
%                 PlotData.ERPbins2(:,:,bin,z,s)= single(squeeze(nanmean(erp(:, :, delay==z & allRT>rtlims3(bin,z,s) & allRT<rtlims3(bin+1,z,s)),3)));
%                 PlotData.ERPrbins2(:,:,bin,z,s)= single(squeeze(nanmean(erpr(:, :, delay==z & allRT>rtlims3(bin,z,s) & allRT<rtlims3(bin+1,z,s)),3)));
%             end
%             
%             
%             numbins2=2;
%             for d=1:5
%                 rtlims2(:,d,z,s) = prctile(allRT(delay==z), linspace(0, 100, numbins2+1)) ;
%                 for bin=1:numbins2
%                     PlotDataTraining.numtrialsBins(bin,d,z,s) = sum(allRT>rtlims2(bin,d,z,s) & allRT<rtlims2(bin+1,d,z,s) &  sessid==d & delay==z & ~Exclude);
%                     PlotDataTraining.RTBins(bin, d,z, s) = nanmean(allRT(allRT>rtlims2(bin,d,z,s) & allRT<rtlims2(bin+1,d,z,s) &  sessid==d & delay==z)); % Get average RT within each RTbin
%                     PlotDataTraining.AccBins(bin,d,z,s) = sum(delay(1,:)==z &  sessid==d & allRT>rtlims2(bin,d,z,s) & allRT<rtlims2(bin+1,d,z,s) & Correct_Error==1)/sum(delay==z &  sessid==d & allRT>rtlims2(bin,d,z,s) & allRT<rtlims2(bin+1,d,z,s));
%                     PlotDataTraining.ERPbins(:,:,bin,d,z,s)= single(squeeze(nanmean(erp(:, :, delay==z &  sessid==d & allRT>rtlims2(bin,d,z,s) & allRT<rtlims2(bin+1,d,z,s)),3)));
%                     PlotDataTraining.ERPrbins(:,:,bin,d,z,s)= single(squeeze(nanmean(erpr(:, :, delay==z & sessid==d & allRT>rtlims2(bin,d,z,s) & allRT<rtlims2(bin+1,d,z,s)),3)));
%                 end
%             end
%             numbins4=4; 
%             for d=1:5
%                 rtlims4(:,d,z,s) = prctile(allRT(delay==z), linspace(0, 100, numbins4+1)) ;
%                 for bin=1:numbins4
%                     PlotDataTraining.numtrialsBins2(bin,d,z,s) = sum(allRT>rtlims4(bin,d,z,s) & allRT<rtlims4(bin+1,d,z,s) &  sessid==d & delay==z & ~Exclude);
%                     PlotDataTraining.RTBins2(bin, d,z, s) = nanmean(allRT(allRT>rtlims4(bin,d,z,s) & allRT<rtlims4(bin+1,d,z,s) &  sessid==d & delay==z)); % Get average RT within each RTbin
%                     PlotDataTraining.AccBins2(bin,d,z,s) = sum(delay(1,:)==z &  sessid==d & allRT>rtlims4(bin,d,z,s) & allRT<rtlims4(bin+1,d,z,s) & Correct_Error==1)/sum(delay==z &  sessid==d & allRT>rtlims4(bin,d,z,s) & allRT<rtlims4(bin+1,d,z,s));
%                     PlotDataTraining.ERPbins2(:,:,bin,d,z,s)= single(squeeze(nanmean(erp(:, :, delay==z &  sessid==d & allRT>rtlims4(bin,d,z,s) & allRT<rtlims4(bin+1,d,z,s)),3)));
%                     PlotDataTraining.ERPrbins2(:,:,bin,d,z,s)= single(squeeze(nanmean(erpr(:, :, delay==z & sessid==d & allRT>rtlims4(bin,d,z,s) & allRT<rtlims4(bin+1,d,z,s)),3)));
%                 end
%             end
%         end
%         PlotData.rtlims=rtlims1; PlotData.rtlims2=rtlims3;PlotDataTraining.rtlims=rtlims2;PlotDataTraining.rtlims2=rtlims4;       
        end
        
 clear erp* erpr* ERPRT ZERP* Exclude AnyArts AnyBlinks allRT allTrig allRTsamp allrespLR delay_Nmin block_Nmin delay delay2 sessid blockid Correct_Error allRTsamp FixLen
        disp(['subject ' allsubj{s} ' HPF ' Files{HPF} ' CSD ' num2str(spline) ' numtrials = ' num2str(numtrialsdisp)]) % display number of valid trials (ie RT greater than 75 but less than 1900
        
        
    end
    
    
    Sessid=single(Sessid); SubID=single(SubID); Block=single(Block); Delay=single(Delay); Targ=single(Targ);
    Resp=single(Resp); RT=single(RT); Acc=single(Acc); RTSamp=single(RTSamp);
    
    %% Save data
    if ismember(spline, 2:7)
        if k==1 
            save([matfolder2 SubFold{HPF} 'SingleTrial_CPP_' Files{HPF} '_CSD' num2str(spline)],'CPP','CPPr','CPPch','artifs','blinks','SubID', 'Sessid', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
            clear CPP CPPr
        elseif k==2 
            save([matfolder2 SubFold{HPF} 'SingleTrial_SPN_' Files{HPF} '_CSD' num2str(spline)],'SPN','SPNr','SPNch','artifs','blinks','SubID', 'Sessid', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
            clear SPN SPNr
        elseif k==3 
            save([matfolder2 SubFold{HPF} 'SingleTrial_LRP_' Files{HPF} '_CSD' num2str(spline)],'LRP','LRPr','LRPch','artifs','blinks','SubID', 'Sessid', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
            clear LRP LRPr
        elseif k==4
            save([matfolder2 SubFold{HPF} 'SingleTrial_FCN_' Files{HPF} '_CSD' num2str(spline) ],'FCN','FCNr','FCNch','artifs','blinks','SubID', 'Sessid', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
            clear FCN FCNr
        end
        save([matfolder2, 'GrAvERP_Sessions' Files{HPF} '_CSD' num2str(spline),'.mat'], 'PlotData', 'PlotDataTraining', 'slock_wind', 'rlock_wind', 'delays');
        clear PlotData PlotDataTraining;
    
    elseif spline==1
        if k==1 
            save([matfolder2 SubFold{HPF} 'SingleTrial_CPP_' Files{HPF} '_NoCSD'],'CPP','CPPr', 'CPPch', 'artifs','blinks','SubID', 'Sessid', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
            clear CPP CPPr
        elseif k==2
            save([matfolder2 SubFold{HPF} 'SingleTrial_SPN_' Files{HPF} '_NoCSD'],'SPN','SPNr', 'SPNch', 'artifs','blinks','SubID', 'Sessid', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
            clear SPN SPNr
        elseif k==3 
            save([matfolder2 SubFold{HPF} 'SingleTrial_LRP_' Files{HPF} '_NoCSD'],'LRP','LRPr', 'LRPch', 'artifs','blinks','SubID', 'Sessid', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
            clear LRP LRPr
        elseif k==4
            save([matfolder2 SubFold{HPF} 'SingleTrial_FCN_' Files{HPF} '_NoCSD'],'FCN','FCNr', 'LRPch', 'artifs','blinks','SubID', 'Sessid', 'Block','Delay', 'Acc','Targ', 'Resp', 'RT', 'slock_wind','slock_winds','rlock_winds');
            clear FCN FCNr
        end
        save([matfolder2, 'GrAvERP_Sessions' Files{HPF} '_NoCSD.mat'], 'PlotData', 'PlotDataTraining','numtrials', 'slock_wind', 'rlock_wind', 'delays');
        clear PlotData PlotDataTraining;
    end
    
    savefig([matfolder2 SubFold{HPF} 'Artifacts_CSD' num2str(spline) '_' num2str(Files{HPF})]);
    
        
end