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


%% Beta 
HPF=3; 
load([folder 'GrAv_TimeFreq_HPF_0_05_CSD4']);
FBeta=ismember(round(Freqs), [10 12 15 17 22 27 30]); %mubeta

BetaChansLR = {[106:110 113:117 122:124] [48:56 61:64]} ;
RTWind=TimePointsR>-150 & TimePointsR<-50;
BLWind=TimePoints>550 & TimePoints<600;

%% TopoPlot & Electrode selection
clear Topo*
for d=1:5
    for z=1:3
        for r=1:2
            TopoRT(:,r,z,d,:) = nanmean(nanmean(PlotTraining1.TFreqr(:,FBeta,RTWind,r,z,d,:),2),3);
            TopoBL(:,r,z,d,:) = nanmean(nanmean(PlotTraining1.TFreq(:,FBeta,BLWind,r,z,d,:),2),3);
        end
    end
end
TopoRTLat= squeeze(TopoRT(:,2,:,:,:)-TopoRT(:,1,:,:,:));

for s=1:length(allsubj)
    for d=1:5
        if s==11 && d==5
            dd=4;
        else dd=d;
        end
    thisdat=nanmean(TopoRT(BetaChansLR{1},2,:,dd,s)-TopoRT(BetaChansLR{1},1,:,dd,s),3); 
    [m idx]=sort(thisdat,'ascend'); 
    choose=sum(isnan(m))+1; 
    Betachans{1}{s}(d,:)= BetaChansLR{1}(idx(choose:choose+3));
    
    thisdat2=nanmean(TopoRT(BetaChansLR{2},1,:,dd,s)-TopoRT(BetaChansLR{2},2,:,dd,s),3);
    [m2 idx2]=sort(thisdat2,'ascend');
    choose2=sum(isnan(m2))+1; % avoid choosing a channel with nans  
    Betachans{2}{s}(d,:)= BetaChansLR{2}(idx2(choose2:choose2+3));
%     zdat=zscore(thisdat+10);
%     outliers=zdat<-3 | zdat>3; 
%     thisdat(outliers)=nan; 
    end
end


TopoDataFFT_rlock=squeeze(PlotTraining1.TFreqr);
TopoDataFFT_slock=squeeze(PlotTraining1.TFreq);
clear Plot*

%% Load single trial data
load(['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\SingleTrialMatrices\LPF_HPF_0_05\SingleTrial_MotorFFT__HPF_0_05_CSD4']);

Delay=Delay';
Delay(2,:)=Delay(2,:)*(1000/512);
RT=RT*(1000/512);
RT2=RT+Delay(2,:);
RT2(RT2<1000)=nan;
rlock_winds = round(-1*fs: .6*fs);
rlock_wind = rlock_winds*(1000/512);
    
FFTMotor(FFTMotor==0)=nan; FFTMotorR(FFTMotorR==0)=nan; 

ZFFT=nan(1,size(FFTMotor,4));

%% Condense Matrices 
evwind=[-50 50]; rwin=[-150 -50];  blwin=BLWind;

for s=1:length(allsubj)
    for d=1:5
        idx=SubID==s & Sessid==d;
        idxch=ismember(MotorChans, Betachans{1}{s}(d,:));
        BL=squeeze(nanmean(nanmean(nanmean(FFTMotor(idxch,FBeta,BLWind,idx),1),2),3))';
        Beta(1,:,idx)=squeeze(nanmean(nanmean(FFTMotor(idxch, FBeta, :, idx))))-repmat(BL, [size(FFTMotor,3) 1]);
        BetaR(1,:,idx)=squeeze(nanmean(nanmean(FFTMotorR(idxch, FBeta, :, idx))))-repmat(BL, [size(FFTMotorR,3) 1]);
        BetaRT(1,idx)=squeeze(nanmean(nanmean(BetaR(1, TimePointsR>rwin(1) & TimePointsR<rwin(2), idx),1),2));
        for z=1:3
            idx2=SubID==s & Sessid==d & Delay(1,:)==z;
            timeidx=TimePoints>delays(z)+evwind(1) & TimePoints<delays(z)+evwind(2);
            BetaEv(1,idx2)=squeeze(nanmean(nanmean(Beta(1, timeidx, idx2),1),2));
        end
            
        idxch2=ismember(MotorChans, Betachans{2}{s}(d,:));
        BL=squeeze(nanmean(nanmean(nanmean(FFTMotor(idxch2,FBeta,BLWind,idx),1),2),3))';
        Beta(2,:,idx)=squeeze(nanmean(nanmean(FFTMotor(idxch2, FBeta, :, idx))))-repmat(BL, [size(FFTMotor,3) 1]);
        BetaR(2,:,idx)=squeeze(nanmean(nanmean(FFTMotorR(idxch2, FBeta, :, idx))))-repmat(BL, [size(FFTMotorR,3) 1]);
        
        BetaRT(2,idx)=squeeze(nanmean(nanmean(BetaR(2, TimePointsR>rwin(1) & TimePointsR<rwin(2), idx),1),2));
        for z=1:3
            idx2=SubID==s & Sessid==d & Delay(1,:)==z;
            timeidx=TimePoints>delays(z)+evwind(1) & TimePoints<delays(z)+evwind(2);
            BetaEv(2,idx2)=squeeze(nanmean(nanmean(Beta(2, timeidx, idx2),1),2));
        end
    end
end

clear FFTMotor FFTMotorR
checknans=isnan(BetaRT(1,:));
AnyBlinks = (blinks.Epoch)>0;
AnyArts = (artifs.Epoch)>0;
Exclude = AnyBlinks>0 | sum(AnyArts,1)>0;
Beta(:,:,Exclude)=nan; BetaR(:,:,Exclude)=nan; BetaEv(:,Exclude)=nan; BetaRT(:,Exclude)=nan; 
save(['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\Analyses\Temporal_Uncertainty_Delays\PlottingData\BetaBL_PlottingData'], 'Exclude', 'Beta*','TimePoints*', 'Freqs', 'FBeta','Topo*', 'SubID', 'Acc', 'Block', 'Delay','Resp', 'RT', 'RT2', 'Sessid', 'Targ')






