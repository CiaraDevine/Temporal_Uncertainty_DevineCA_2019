%% %%%%%%%%%%%%%%%% Separate Delays %%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

labels = {'No HPF' 'HPF 0.05' 'HPF 0.15' 'HPF 0.25' 'HPF 0.35' 'HPF 0.45'};
Files = {['_NoHPF'] ['_HPF_0_05'] ['_HPF_0_15'] ['_HPF_0_25'] ['_HPF_0_35'] ['_HPF_0_45']};
folder = ['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\SingleTrialMatrices\'];


allsubj =  { 'CJG'  'CR'  'AR' 'KD' 'TB'  'AOB'  'PM'  'CMG'  'SR'  'AB' 'GK' 'ED' 'SC' 'CB' 'CE' 'SH' 'RMC' 'JH' 'JR' 'ROC' 'KM' 'ROS'}; % In order of when they were tested %JR and KD excluded
sessions = {[1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:4] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5] [1:5]  [1:5]};

fs=512;

colours_sessions = {[1, .65, .67] [.95, .35, .40] [1, 0, .06] [.5, 0, 0] [0, 0, 0]}; % Red theme colours
colours_blocks ={[.46, 1, .58] [.28, 1, .4] [.28, 1, .88] [.14, .97, 1] [0, .86, .94] [.08, .63, 1] [.02, .4, 1] [.28, 0, .8] [.01, 0, .4] [0, 0, 0]}; %green/blue themed colours
chanlocs = readlocs('cap128.loc');
delays=[800 1200 1600];
scrsz = get(0,'ScreenSize');

subs=[2:5 7:18 20:22];


%% SSVEP 
HPF=3; 
load([folder 'GrAv_TimeFreq_HPF_0_05_CSD4']);
FBeta=ismember(round(Freqs), [12 15 17 22 27 30]); %mubeta

F20=ismember(round(Freqs,1), [20]); 
F20_Adj = ismember(round(Freqs,1), [17.5 22.5 27.5]);
F25=ismember(round(Freqs,1), [25]); 
F25_Adj = ismember(round(Freqs,1), [17.5 22.5 27.5]);

RTWind=TimePointsR>-50 & TimePointsR<50;
BLWind=TimePoints>500 & TimePoints<600;

%% TopoPlot & Electrode selection
for z=1:3
    for r=1:2
        Topo20RT(:,r,z,:) = nanmean(nanmean(PlotData.TFreqr(:,[F20],RTWind,r,z,:),2),3);%./nanmeannanmean(PlotData.TFreqr(:,[F20_Adj],RTWind,r,z,:),2),3);
        Topo25RT(:,r,z,:) = nanmean(nanmean(PlotData.TFreqr(:,[F25],RTWind,r,z,:),2),3);%./nanmeannanmean(PlotData.TFreqr(:,[F25_Adj],RTWind,r,z,:),2),3);
        if r==1, TopoDiff(:,r,z,:) = Topo20RT(:,r,z,:)-Topo25RT(:,r,z,:);
        else TopoDiff(:,r,z,:) = Topo25RT(:,r,z,:)-Topo20RT(:,r,z,:); end
    end
end

ylims=[-.18 .18];
%for s=1:15
f1=figure;
title('Topo RT', 'fontsize', 22, 'fontweight', 'bold');
h=topoplot(nanmean(nanmean(nanmean(TopoDiff(:,:,:,:),2),3),4), chanlocs,'plotchans',[1:128], ...
    'maplimits', ylims,'electrodes','off', 'whitebk', 'off', 'headrad',.65, 'plotrad', .66, 'emarker2', {[],'o','k',10,1}, 'hcolor', 'k', 'conv', 'on');
%set(h, 'position', [.2 .2 1.5 1.5], 'outerposition', [0 0 1 1]);
cb= colorbar; set(cb, 'location', 'southoutside', 'fontsize', 22, 'position', [.73 .12 .26 .06]); title(cb, ' (uVm ^2)');
%set(f1,'position',scrsz);
%end
OccChans = [ 5 14:24 27:32];
for s=1:length(allsubj)
    thisdat=nanmean(nanmean(TopoDiff(OccChans,:,:,s),3),2);
    zdat=zscore(thisdat+1);
    outliers=zdat<-3 | zdat>3; 
    thisdat(outliers)=nan; 
    [m idx]=sort(thisdat,'descend');
    choose=sum(isnan(m))+1; % avoid choosing a channel with nans 
    SSVEPchans{s}= OccChans(idx(choose:choose+3));
end

% TopoDataFFT_rlock=squeeze(nanmean(PlotData.TFreqr,4));
% TopoDataFFT_slock=squeeze(nanmean(PlotData.TFreq,4));

 clear Plot*

 
%% Load single trial data
load([folder 'LPF_HPF_0_05\SingleTrial_OccFFT__HPF_0_05_CSD4']);
Delay=Delay';
Delay(2,:)=Delay(2,:)*(1000/512);
RT=RT*(1000/512);
RT2=RT+Delay(2,:);
RT2(RT2<1000)=nan;
rlock_winds = round(-1*fs: .6*fs);
rlock_wind = rlock_winds*(1000/512);
TimePoints=TimePoints(1:size(FFTOcc,3));
TimePointsR=TimePointsR(1:size(FFTOccR,3));

%% Exclude Artifact/Blink Trials
AnyBlinks = (blinks.Epoch)>0;
AnyArts = (artifs.Epoch)>0;
Exclude = AnyBlinks>0 | sum(AnyArts,1)>0;
clear blinks artifs AnyArts AnyBlinks

SSVEP20Hz=nan(size(FFTOcc,3),size(FFTOcc,4)); 
SSVEP25Hz=nan(size(FFTOcc,3),size(FFTOcc,4)); 
SSVEP25HzAdj=nan(size(FFTOcc,3),size(FFTOcc,4)); 
SSVEP20HzAdj=nan(size(FFTOcc,3),size(FFTOcc,4)); 
SSVEPr20Hz=nan(size(FFTOccR,3),size(FFTOcc,4)); 
SSVEPr25Hz=nan(size(FFTOccR,3),size(FFTOcc,4)); 
SSVEPr25HzAdj=nan(size(FFTOccR,3),size(FFTOcc,4)); 
SSVEPr20HzAdj=nan(size(FFTOccR,3),size(FFTOcc,4));

ZFFT=nan(1,size(FFTOcc,4));
FFTOcc(FFTOcc==0)=nan; FFTOccR(FFTOccR==0)=nan; 
FFTOccR(:,:,:,Exclude)=nan; FFTOcc(:,:,:,Exclude)=nan; 

for s=1:length(allsubj)
    if s~=19
        for d=sessions{s}
            idx=SubID==s & Sessid==d;
            idxch=ismember(OccChans, SSVEPchans{s});
            SSVEP20Hz(:,idx)=squeeze(nanmean(nanmean(FFTOcc(idxch, F20, :,idx),1),2));
            SSVEP20HzAdj(:,idx)=squeeze(nanmean(nanmean(FFTOcc(idxch, F20_Adj, :,idx),1),2));
            SSVEP25Hz(:,idx)=squeeze(nanmean(nanmean(FFTOcc(idxch, F25, :,idx),1),2));
            SSVEP25HzAdj(:,idx)=squeeze(nanmean(nanmean(FFTOcc(idxch, F25_Adj, :,idx),1),2));
            
            SSVEPr20Hz(:,idx)=squeeze(nanmean(nanmean(FFTOccR(idxch, F20, :,idx),1),2));
            SSVEPr20HzAdj(:,idx)=squeeze(nanmean(nanmean(FFTOccR(idxch, F20_Adj, :,idx),1),2));
            SSVEPr25Hz(:,idx)=squeeze(nanmean(nanmean(FFTOccR(idxch, F25, :,idx),1),2));
            SSVEPr25HzAdj(:,idx)=squeeze(nanmean(nanmean(FFTOccR(idxch, F25_Adj, :,idx),1),2));
        end
    end
end

SSVEPSNR_20Hz = SSVEP20Hz./SSVEP20HzAdj; 
SSVEPSNRr_20Hz = SSVEPr20Hz./SSVEPr20HzAdj; 
SSVEPSNR_25Hz = SSVEP25Hz./SSVEP25HzAdj; 
SSVEPSNRr_25Hz = SSVEPr25Hz./SSVEPr25HzAdj;

save(['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\Analyses\Temporal_Uncertainty_Delays\PlottingData\SSVEP_PlottingData'], 'SSVEPSNR*', 'TimePoints', 'TimePointsR',  'Topo*', 'delays','SubID', 'Acc', 'Block', 'Delay','Resp', 'RT', 'RT2', 'Sessid', 'Targ', 'subs', 'colours*')

%%
 windR=TimePointsR>-300 & TimePointsR<-200; 

for s=1:length(allsubj)
    for z=1:3
        idxLt= Targ==8 & Delay(1,:)==z & SubID==s; idxRt= Targ==9 & Delay(1,:)==z & SubID==s & RT>500; 
        idxLr= Resp==12 & Delay(1,:)==z & SubID==s; idxRr= Resp==13 & Delay(1,:)==z & SubID==s & RT>500; 
        SSVEP20Lt(:,z,s)=nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20, :,idxLt),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,idxLt),1),2),4);
        SSVEP20Rt(:,z,s)=nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20, :,idxRt),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,idxRt),1),2),4);
        SSVEP25Lt(:,z,s)=nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25, :,idxLt),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,idxLt),1),2),4);
        SSVEP25Rt(:,z,s)=nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25, :,idxRt),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,idxRt),1),2),4);
        SSVEPr20Lt(:,z,s)=nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20, :,idxLt),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,idxLt),1),2),4);
        SSVEPr20Rt(:,z,s)=nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20, :,idxRt),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,idxRt),1),2),4);
        SSVEPr25Lt(:,z,s)=nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25, :,idxLt),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,idxLt),1),2),4);
        SSVEPr25Rt(:,z,s)=nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25, :,idxRt),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,idxRt),1),2),4);
           
        numbins=8;
        RTBinsLt(:,z,s)=prctile(RT2(idxLt & RT<1800), linspace(0,100,numbins+1)); RTBinsRt(:,z,s)=prctile(RT2(idxRt & RT<1800), linspace(0,100,numbins+1));
        RTBinsLr(:,z,s)=prctile(RT2(idxLr & RT<1800), linspace(0,100,numbins+1)); RTBinsRr(:,z,s)=prctile(RT2(idxRr & RT<1800), linspace(0,100,numbins+1));
        for b=1:numbins
            
            idxRTLt = RT2>RTBinsLt(b,z,s) & RT2<RTBinsLt(b+1,z,s) & idxLt; numtrialsLt(b,z,s)=sum(idxRTLt);
            idxRTRt = RT2>RTBinsRt(b,z,s) & RT2<RTBinsRt(b+1,z,s) & idxRt; numtrialsRt(b,z,s)=sum(idxRTRt);
            idxRTLr = RT2>RTBinsLr(b,z,s) & RT2<RTBinsLr(b+1,z,s) & idxLr; numtrialsLr(b,z,s)=sum(idxRTLr);
            idxRTRr = RT2>RTBinsRr(b,z,s) & RT2<RTBinsRr(b+1,z,s) & idxRr; numtrialsRr(b,z,s)=sum(idxRTRr);
            RTBins(b,z,s)=nanmean(RT2(idxRTLt | idxRTRt));
            
            ZFFT(idxRTLt)=squeeze(zscore(nanmean(nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20, windR, idxRTLt),1),2),3),4)./nanmean(nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20_Adj, windR, idxRTLt),1),2),3),4)));
            ZFFT(idxRTRt)=squeeze(zscore(nanmean(nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25, windR, idxRTRt),1),2),3),4)./nanmean(nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25_Adj, windR, idxRTRt),1),2),3),4)));
            Outliers= ZFFT<-3 & ZFFT>3;
            
            SSVEP20Lt2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20, :,  idxRTLt & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,  idxRTLt & ~Outliers),1),2),4);
            SSVEP20Rt2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20, :,  idxRTRt & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,  idxRTRt & ~Outliers),1),2),4);
            SSVEP25Lt2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25, :,  idxRTLt & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,  idxRTLt & ~Outliers),1),2),4);
            SSVEP25Rt2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25, :,  idxRTRt & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,  idxRTRt & ~Outliers),1),2),4);
            
            SSVEPr20Lt2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20, :,  idxRTLt & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,  idxRTLt & ~Outliers),1),2),4);
            SSVEPr20Rt2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20, :,  idxRTRt & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,  idxRTRt & ~Outliers),1),2),4);
            SSVEPr25Lt2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25, :,  idxRTLt & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,  idxRTLt & ~Outliers),1),2),4);
            SSVEPr25Rt2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25, :,  idxRTRt & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,  idxRTRt & ~Outliers),1),2),4);
            
            SSVEP20Lr2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20, :,  idxRTLr & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,  idxRTLr & ~Outliers),1),2),4);
            SSVEP20Rr2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20, :,  idxRTRr & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,  idxRTRr & ~Outliers),1),2),4);
            SSVEP25Lr2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25, :,  idxRTLr & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,  idxRTLr & ~Outliers),1),2),4);
            SSVEP25Rr2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25, :,  idxRTRr & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOcc(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,  idxRTRr & ~Outliers),1),2),4);
                       
            SSVEPr20Lr2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20, :,  idxRTLr & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,  idxRTLr & ~Outliers),1),2),4);
            SSVEPr20Rr2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20, :,  idxRTRr & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F20_Adj, :,  idxRTRr & ~Outliers),1),2),4);
            SSVEPr25Lr2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25, :,  idxRTLr & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,  idxRTLr & ~Outliers),1),2),4);
            SSVEPr25Rr2(:,b,z,s) = nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25, :,  idxRTRr & ~Outliers),1),2),4)./nanmean(nanmean(nanmean(FFTOccR(ismember(OccChans, SSVEPchans{s}), F25_Adj, :,  idxRTRr & ~Outliers),1),2),4);
            

        end
    end
end


%% Stats 
TargSSVEPdat=(SSVEP20Lt(:,:,subs)+SSVEP25Rt(:,:,subs))/2;
NonTargSSVEPdat=(SSVEP20Rt(:,:,subs)+SSVEP25Lt(:,:,subs))/2; 
for z=1:3
Inputdata(1,z,:)=nanmean(TargSSVEPdat(TimePoints>delays(z)+400 & TimePoints<delays(z)+600, z,:)); 
Inputdata(2,z,:)=nanmean(NonTargSSVEPdat(TimePoints>delays(z)+400 & TimePoints<delays(z)+600, z, :));
end
ANOVA.SSVEPPostEv=RMAOV(Inputdata);

TargSSVEPrdat=(SSVEPr20Lt(:,:,subs)+SSVEPr25Rt(:,:,subs))/2;
NonTargSSVEPrdat=(SSVEPr20Rt(:,:,subs)+SSVEPr25Lt(:,:,subs))/2; 
Inputdata(1,z,:)=nanmean(TargSSVEPrdat(TimePointsR>-250 & TimePointsR<-200, z,:)); 
Inputdata(2,z,:)=nanmean(NonTargSSVEPrdat(TimePointsR>-250 & TimePointsR<-200,z, :));
ANOVA.SSVEPPreRT=RMAOV(Inputdata);

%%     
 figure; 
ylims=[-.4 .8]; hold on;
for z=1:3
    subplot(1,3,z); hold on; 
    for b=1:numbins
    TargSSVEP=((nanmean(SSVEP20Lr2(:,b,z,subs),4)+nanmean(SSVEP25Rr2(:,b,z,subs),4))/2);
    NonTargSSVEP=((nanmean(SSVEP20Rr2(:,b,z,subs),4)+nanmean(SSVEP25Lr2(:,b,z,subs),4))/2);  
    Xaxisidx=TimePoints>500 & TimePoints<(delays(z));
    Plc(z)=plot(TimePoints(Xaxisidx), cumsum(TargSSVEP(Xaxisidx)-NonTargSSVEP(Xaxisidx)),'linewidth', 4, 'linestyle', '-',  'color', colours_blocks{b+1});
    %Pli(z)=plot(TimePoints(Xaxisidx), cumsum(NonTargSSVEP(Xaxisidx)),'linewidth', 4, 'linestyle', '--',  'color', colours_blocks{b+1});
    xlim([0 3400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
    line([delays(z) delays(z)],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 20);
    end
end
 figure; 
ylims=[1.2 1.8]; hold on;
for z=1:3
    subplot(1,3,z); hold on; 
    for b=1:numbins
    TargSSVEPr=((nanmean(SSVEPr20Lt2(:,b,z,subs),4)+nanmean(SSVEPr25Rt2(:,b,z,subs),4))/2);
    NonTargSSVEPr=((nanmean(SSVEPr20Rt2(:,b,z,subs),4)+nanmean(SSVEPr25Lt2(:,b,z,subs),4))/2);    
    Plc(z)=plot(TimePointsR, TargSSVEPr,'linewidth', 4, 'linestyle', '-',  'color', colours_blocks{b+1});
    Pli(z)=plot(TimePointsR, NonTargSSVEPr,'linewidth', 4, 'linestyle', '--',  'color', colours_blocks{b+1});
    xlim([-800 400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
    line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 20, 'ycolor', 'w', 'yticklabels', [], 'ytick', [],'xtick', [-800 -400 0 400], 'xticklabels', {'-800' '-400' '0' '400'});
    end
end



%%     
 figure; 
ylims=[-3 8]; hold on;
for z=1:3
    subplot(1,3,z); hold on; 
    for b=1:numbins
    TargSSVEP=((nanmean(SSVEP20Lt2(:,b,z,subs),4)+nanmean(SSVEP25Rt2(:,b,z,subs),4))/2);
    NonTargSSVEP=((nanmean(SSVEP20Rt2(:,b,z,subs),4)+nanmean(SSVEP25Lt2(:,b,z,subs),4))/2); 
    SSVEPDiff=TargSSVEP-NonTargSSVEP;
    Plc(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), cumsum(SSVEPDiff(TimePoints<(delays(z)+1800))),'linewidth', 4, 'linestyle', '-',  'color', colours_blocks{b+1});
    xlim([0 3400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
    line([delays(z) delays(z)],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 20);
    end
end
 figure; 
ylims=[-2 5.8]; hold on;
for z=1:3
    subplot(1,3,z); hold on; 
    for b=1:numbins
    TargSSVEPr=((nanmean(SSVEPr20Lt2(:,b,z,subs),4)+nanmean(SSVEPr25Rt2(:,b,z,subs),4))/2);
    NonTargSSVEPr=((nanmean(SSVEPr20Rt2(:,b,z,subs),4)+nanmean(SSVEPr25Lt2(:,b,z,subs),4))/2);    
   SSVEPDiffr=TargSSVEPr-NonTargSSVEPr;
    Plc(z)=plot(TimePointsR, cumsum(SSVEPDiffr),'linewidth', 4, 'linestyle', '-',  'color', colours_blocks{b+1});
    xlim([-800 400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
    line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 20, 'ycolor', 'w', 'yticklabels', [], 'ytick', [],'xtick', [-800 -400 0 400], 'xticklabels', {'-800' '-400' '0' '400'});
    end
end

%% Plot SSVEP
 
for z=1:3
    wind(:,z)=TimePoints>delays(z)-300 & TimePoints<delays(z)-100; 
end
 windR=TimePointsR>-600 & TimePointsR<-200; 
 
figure; 
ylims=[1 1.6];
sp(3)= subplot(2,5,3); hold on;
for z=1:3
    TargSSVEPr=((nanmean(SSVEPr20Lt(:,z,subs),3)+nanmean(SSVEPr25Rt(:,z,subs),3))/2);
    NonTargSSVEPr=((nanmean(SSVEPr20Rt(:,z,subs),3)+nanmean(SSVEPr25Lt(:,z,subs),3))/2);    
    Plc(z)=plot(TimePointsR, TargSSVEPr,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    Pli(z)=plot(TimePointsR, NonTargSSVEPr ,'linewidth', 4, 'linestyle', '--',  'color', colours_sessions{z+2});
    xlim([-810 400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
    line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 20, 'ycolor', 'w', 'yticklabels', [], 'ytick', [],'xtick', [-800 -400 0 400], 'xticklabels', {'-800' '-400' '0' '400'});
end


sp(2)=subplot(2,5,2); hold on;
%title(['CPP ' labels{HPF} ' CSD SPline ' num2str(spline)], 'fontsize', 20, 'fontweight', 'bold');
for z=1:3
    TargSSVEP=((nanmean(SSVEP20Lt(TimePoints<(delays(z)+1800),z,subs),3)+nanmean(SSVEP25Rt(TimePoints<(delays(z)+1800),z,subs),3))/2);
    NonTargSSVEP=((nanmean(SSVEP20Rt(TimePoints<(delays(z)+1800),z,subs),3)+nanmean(SSVEP25Lt(TimePoints<(delays(z)+1800),z,subs),3))/2);    
    Plc(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), TargSSVEP ,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    Plc1(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), NonTargSSVEP ,'linewidth', 4, 'linestyle', '--',  'color', colours_sessions{z+2});
    %Amp=nanmean((slock_wind{z}>delays(z) & slock_wind{z}<delays(z)+25));
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    
end
xlabel('Time (ms)', 'fontsize', 18);
xlim([500 3600]); ylim(ylims);
line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
set(gca, 'fontsize', 20, 'xtick', [800: 400 : 3600], 'xticklabels', {'800' '1200' '1600' '2000' '2400' '2800' '3200' ''});
%ylabel('CPP Amplitude uV', 'fontsize', 26);
ylabel('SSVEP SNR ', 'fontsize', 26);

% lg=legend([Plc1(3) Pli(3)], 'Contra', 'Ipsi');
% set(lg, 'fontsize', 26, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast', 'box', 'off');
% % h=axes; set(h,'position', [0 .5 .49 .49], 'visible', 'off');label=text(.02, .97, 'A'); set(label,'fontsize', 34, 'fontweight', 'bold');
% % h1=axes; set(h1,'position', [.51 .5 .49 .49], 'visible', 'off');label=text(.3, .97, 'B'); set(label,'fontsize', 34, 'fontweight', 'bold');


ylims=[-.51 .51];
sp(1)=subplot(2,5,1); topoplot(nanmean(nanmean(nanmean(TopoDiff,2),3),4), chanlocs,'plotchans',[1:128],'maplimits', ylims, ...
    'electrodes','off', 'whitebk', 'off', 'headrad',.65, 'plotrad', .66, 'emarker2', {[],'o','k',10,1}, 'hcolor', 'k', 'conv', 'on');
%set(h, 'position', [.2 .2 1.5 1.5], 'outerposition', [0 0 1 1]);
%cb= colorbar; set(cb, 'location', 'eastoutside', 'fontsize', 22, 'position', [.95 .63 .015 .25]); %title(cb, ' (uVm ^2)');

ylims=[1.1 1.8]; 
ylims=[-.2 .6];
sp(6)=subplot(2,5,6); hold on; box off; 
for z=1:3
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),abs(nanmean(BetaRContraBins(:,z,subs),3)-nanmean(BetaRIpsiBins(:,z,subs),3)),'linestyle', '-', 'linewidth',3, 'color', colours_sessions{z+2});
    TargSSVEPr=((nanmean(nanmean(SSVEPr20Lt2(windR,:,z,subs),1),4)+nanmean(nanmean(SSVEPr25Rt2(windR,:,z,subs),1),4))/2);
    NonTargSSVEPr=((nanmean(nanmean(SSVEPr20Rt2(windR,:,z,subs),1),4)+nanmean(nanmean(SSVEPr25Lt2(windR,:,z,subs),1),4))/2);    
    plot(nanmean(RTBins(:,z,subs),3),  TargSSVEPr,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    plot(nanmean(RTBins(:,z,subs),3), NonTargSSVEPr,'linewidth', 4, 'linestyle', '--',  'color', colours_sessions{z+2});
    plot(nanmean(RTBins(:,z,subs),3), TargSSVEPr-NonTargSSVEPr,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [500 800 1200 1600 2000 2400 2800 3200], 'xticklabel', {'0' ' ..800' '1200' '1600' '2000' '2400' '2800' '3600'});
ylabel('SSVEP SNR @ RT','fontsize',22); xlabel('RT (ms)', 'fontsize',22);
% lg=legend([pl(3) pl2(3)],'At Response', 'Pre Contrast Change');
% set(lg, 'fontsize', 24, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast');

sp(7)=subplot(2,5,7); hold on; box off; 
for z=1:3
    TargSSVEP=((nanmean(nanmean(SSVEP20Lt2(wind(:,z),:,z,subs),1),4)+nanmean(nanmean(SSVEP25Rt2(wind(:,z),:,z,subs),1),4))/2);
    NonTargSSVEP=((nanmean(nanmean(SSVEP20Rt2(wind(:,z),:,z,subs),1),4)+nanmean(nanmean(SSVEP25Lt2(wind(:,z),:,z,subs),1),4))/2);    
    Plc(z)=plot(nanmean(RTBins(:,z,subs),3),  TargSSVEP,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    Pli(z)=plot(nanmean(RTBins(:,z,subs),3), NonTargSSVEP ,'linewidth', 4, 'linestyle', '--',  'color', colours_sessions{z+2});
    plot(nanmean(RTBins(:,z,subs),3), TargSSVEP-NonTargSSVEP,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
     l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [], 'xticklabel', {});
ylabel('SSVEP SNR @ Contrast Change','fontsize',22); %xlabel('RT (ms)', 'fontsize',24);

sp(8)=subplot(2,5,8); hold on; box off; 
for z=1:3
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),abs(nanmean(BetaRContraBins(:,z,subs),3)-nanmean(BetaRIpsiBins(:,z,subs),3)),'linestyle', '-', 'linewidth',3, 'color', colours_sessions{z+2});
    TargSSVEPr=((nanmean(nanmean(SSVEPr20Lr2(windR,:,z,subs),1),4)+nanmean(nanmean(SSVEPr25Rr2(windR,:,z,subs),1),4))/2);
    NonTargSSVEPr=((nanmean(nanmean(SSVEPr20Rr2(windR,:,z,subs),1),4)+nanmean(nanmean(SSVEPr25Lr2(windR,:,z,subs),1),4))/2);    
    Plc(z)=plot(nanmean(RTBins(:,z,subs),3),  TargSSVEPr,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    Pli(z)=plot(nanmean(RTBins(:,z,subs),3), NonTargSSVEPr ,'linewidth', 4, 'linestyle', '--',  'color', colours_sessions{z+2});
    plot(nanmean(RTBins(:,z,subs),3), TargSSVEPr-NonTargSSVEPr,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});    
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [500 800 1200 1600 2000 2400 2800 3200], 'xticklabel', {'0' ' ..800' '1200' '1600' '2000' '2400' '2800' '3600'});
ylabel('SSVEP SNR @ RT','fontsize',22); xlabel('RT (ms)', 'fontsize',22);
% lg=legend([pl(3) pl2(3)],'At Response', 'Pre Contrast Change');
% set(lg, 'fontsize', 24, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast');

sp(9)=subplot(2,5,9); hold on; box off; 
for z=1:3
    TargSSVEP=((nanmean(nanmean(SSVEP20Lr2(wind(:,z),:,z,subs),1),4)+nanmean(nanmean(SSVEP25Rr2(wind(:,z),:,z,subs),1),4))/2);
    NonTargSSVEP=((nanmean(nanmean(SSVEP20Rr2(wind(:,z),:,z,subs),1),4)+nanmean(nanmean(SSVEP25Lr2(wind(:,z),:,z,subs),1),4))/2);    
    Plc(z)=plot(nanmean(RTBins(:,z,subs),3),  TargSSVEP,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    Pli(z)=plot(nanmean(RTBins(:,z,subs),3),NonTargSSVEP,'linewidth', 4, 'linestyle', '--',  'color', colours_sessions{z+2});
    plot(nanmean(RTBins(:,z,subs),3), TargSSVEP-NonTargSSVEP,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [], 'xticklabel', {});
ylabel('SSVEP SNR @ Contrast Change','fontsize',22); %xlabel('RT (ms)', 'fontsize',24);

set(sp(1), 'position', [-.05 .7 .25 .25]);
set(sp(2), 'position', [.2 .59 .3 .4]);
set(sp(3), 'position', [.52 .59 .12 .4]);
set(sp(6), 'position', [.08 .08 .3 .2]);
set(sp(7), 'position', [.08 .31 .3 .2]);
set(sp(8), 'position', [.48 .08 .3 .2]);
set(sp(9), 'position', [.48 .31 .3 .2]);

h=axes; set(h,'position', [0 .01 .99 .99], 'visible', 'off');% Create a new axis that is invisible. Otherwise all text will be positioned relevant to your most recent lot which you might not want.
text('String', 'Sensory Representation','fontsize', 30, 'fontweight', 'bold', 'position',[.01, .94]);
text('String', 'F','fontsize', 34, 'fontweight', 'bold', 'position',[.01, .78]);
text('String', 'G','fontsize', 34, 'fontweight', 'bold', 'position',[.16, .97]);
text('String', 'H', 'fontsize', 34, 'fontweight', 'bold','position',[.51, .97]);
text('String', 'I','fontsize', 34, 'fontweight', 'bold','position',[.68, .97]);
text('String','J','fontsize', 34, 'fontweight', 'bold','position',[.68, .72]);



%% Thesis Plot 

 
for z=1:3
    wind(:,z)=TimePoints>delays(z)-300 & TimePoints<delays(z)-100; 
end
 windR=TimePointsR>-600 & TimePointsR<-200; 
 
 
 
figure; 

ylims=[-.61 .51];
sp(1)=subplot(2,3,1); topoplot(nanmean(nanmean(nanmean(TopoDiff,2),3),4), chanlocs,'plotchans',[1:128],'maplimits', ylims, ...
    'electrodes','off', 'whitebk', 'off', 'headrad',.65, 'plotrad', .66, 'emarker2', {[],'o','k',10,1}, 'hcolor', 'k', 'conv', 'on');
%set(h, 'position', [.2 .2 1.5 1.5], 'outerposition', [0 0 1 1]);
%cb= colorbar; set(cb, 'location', 'eastoutside', 'fontsize', 22, 'position', [.95 .63 .015 .25]); %title(cb, ' (uVm ^2)');

ylims=[1.11 1.49];
sp(3)= subplot(2,3,3); hold on;
for z=1:3
    TargSSVEPr=((nanmean(SSVEPr20Lt(:,z,subs),3)+nanmean(SSVEPr25Rt(:,z,subs),3))/2);
    NonTargSSVEPr=((nanmean(SSVEPr20Rt(:,z,subs),3)+nanmean(SSVEPr25Lt(:,z,subs),3))/2);    
    Plc(z)=plot(TimePointsR, TargSSVEPr,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    Pli(z)=plot(TimePointsR, NonTargSSVEPr ,'linewidth', 4, 'linestyle', '--',  'color', colours_sessions{z+2});
    xlim([-810 400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
    line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 20, 'ycolor', 'w', 'yticklabels', [], 'ytick', [],'xtick', [-800 -400 0 400], 'xticklabels', {'' '' '' ''});
end


sp(2)=subplot(2,3,2); hold on;
%title(['CPP ' labels{HPF} ' CSD SPline ' num2str(spline)], 'fontsize', 20, 'fontweight', 'bold');
for z=1:3
    TargSSVEP=((nanmean(SSVEP20Lt(TimePoints<(delays(z)+1800),z,subs),3)+nanmean(SSVEP25Rt(TimePoints<(delays(z)+1800),z,subs),3))/2);
    NonTargSSVEP=((nanmean(SSVEP20Rt(TimePoints<(delays(z)+1800),z,subs),3)+nanmean(SSVEP25Lt(TimePoints<(delays(z)+1800),z,subs),3))/2);    
    Plc(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), TargSSVEP ,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    Plc1(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), NonTargSSVEP ,'linewidth', 4, 'linestyle', '--',  'color', colours_sessions{z+2});
    %Amp=nanmean((slock_wind{z}>delays(z) & slock_wind{z}<delays(z)+25));
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    
end
%xlabel('Time (ms)', 'fontsize', 18);
xlim([500 3600]); ylim(ylims);
line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
set(gca, 'fontsize', 20, 'xtick', [800: 400 : 3600], 'xticklabels', {});
%ylabel('CPP Amplitude uV', 'fontsize', 26);
ylabel('SSVEP SNR ', 'fontsize', 26);


ylims=[-.049 .29];
sp(6)= subplot(2,3,6); hold on;
for z=1:3
    TargSSVEPr=((nanmean(SSVEPr20Lt(:,z,subs),3)+nanmean(SSVEPr25Rt(:,z,subs),3))/2);
    NonTargSSVEPr=((nanmean(SSVEPr20Rt(:,z,subs),3)+nanmean(SSVEPr25Lt(:,z,subs),3))/2);    
    Plc(z)=plot(TimePointsR, TargSSVEPr-NonTargSSVEPr,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    xlim([-810 400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
    line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 20, 'ycolor', 'w', 'yticklabels', [], 'ytick', [],'xtick', [-800 -400 0 400], 'xticklabels', {'-800' '-400' '0' '400'});
end


sp(5)=subplot(2,3,5); hold on;
%title(['CPP ' labels{HPF} ' CSD SPline ' num2str(spline)], 'fontsize', 20, 'fontweight', 'bold');
for z=1:3
    TargSSVEP=((nanmean(SSVEP20Lt(TimePoints<(delays(z)+1800),z,subs),3)+nanmean(SSVEP25Rt(TimePoints<(delays(z)+1800),z,subs),3))/2);
    NonTargSSVEP=((nanmean(SSVEP20Rt(TimePoints<(delays(z)+1800),z,subs),3)+nanmean(SSVEP25Lt(TimePoints<(delays(z)+1800),z,subs),3))/2);    
    Plc(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), TargSSVEP-NonTargSSVEP ,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    %Amp=nanmean((slock_wind{z}>delays(z) & slock_wind{z}<delays(z)+25));
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    
end
xlabel('Time (ms)', 'fontsize', 18);
xlim([500 3600]); ylim(ylims);
line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
set(gca, 'fontsize', 20, 'xtick', [800: 400 : 3600], 'xticklabels', {'800' '1200' '1600' '2000' '2400' '2800' '3200' ''});
%ylabel('CPP Amplitude uV', 'fontsize', 26);
ylabel('dSSVEP ', 'fontsize', 26);



set(sp(1), 'position', [-.05 .5 .25 .25]);
set(sp(2), 'position', [.25 .58 .4 .4]);
set(sp(3), 'position', [.68 .58 .22 .4]);
set(sp(5), 'position', [.25 .08 .4 .4]);
set(sp(6), 'position', [.68 .08 .22 .4]);

