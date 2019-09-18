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
%% Make 6hz filter for visualising CPP
[B,A]=butter(4,(6*2)/fs); % creating filter for ERP waveforms


%% Select Electrodes 

load([folder 'GrAvERP_Sessions_HPF_0_05_CSD4']);

rlock_winds = round(-1*fs: .6*fs);
rlock_wind = rlock_winds*(1000/512);

%% Choose Baseline Correction Window
BL1=slock_wind{3}>-100 & slock_wind{3}<0;
BL2=slock_wind{3}>500 & slock_wind{3}<550; 

BL=BL2; 

PlotDataTraining.ERPPreBL=PlotDataTraining.ERP;
PlotDataTraining.ERPrPreBL=PlotDataTraining.ERPr;

Dim=size(PlotDataTraining.ERPPreBL,2);

PlotDataTraining.ERP=PlotDataTraining.ERPPreBL-repmat(nanmean(PlotDataTraining.ERPPreBL(:,BL,:,:,:,:),2),[1 Dim 1 1 1 1]);
PlotDataTraining.ERPr=PlotDataTraining.ERPrPreBL-repmat(nanmean(PlotDataTraining.ERPPreBL(:,BL,:,:,:,:),2),[1 length(rlock_wind) 1 1 1 1]);

%% TopoPlots

for d=1:5
    TopoRT(:,d,:) = nanmean(nanmean(nanmean(PlotDataTraining.ERPr(:,rlock_wind>-150 & rlock_wind<-50,d,:,1:3,:),2),4),5);
    TopoDiff(:,d,:) = nanmean(nanmean(nanmean(PlotDataTraining.ERPr(:,rlock_wind>-50 & rlock_wind<50,d,:,3,:)-PlotDataTraining.ERPr(:,rlock_wind>-50 & rlock_wind<50,d,:,1,:),2),4),5);
    
    Plotop(:,:,d,:) = nanmean(PlotDataTraining.ERP(:,:,d,:,3,:),4);
    Plotopr(:,:,d,:) = nanmean(nanmean(PlotDataTraining.ERPr(:,:,d,:,:,:),4),5);
end

%ylims=[min(nanmean(nanmean(TopoRT,2),3))-3 max(nanmean(nanmean(TopoRT,2),3))+3]; %

ylims=[-10 10];
f1=figure;
% title(['Topo RT Spline ' num2str(spline) ' HPF ' labels{HPF} '_BL' num2str(BL)], 'fontsize', 22, 'fontweight', 'bold');
%ylims=[min(nanmean(nanmean(TopoEv,2),3))-1 max(nanmean(nanmean(TopoEv,2),3))+1]; ylims= [-20 20];
h=topoplot(nanmean(nanmean(TopoRT(:,:,subs),2),3), chanlocs,'plotchans',[1:128],'maplimits', ylims, ...
'electrodes','labels', 'whitebk', 'off', 'headrad',.65, 'plotrad', .71, 'emarker2', {[3 4],'o','k',10,1}, 'hcolor', 'k', 'conv', 'on');
%set(h, 'position', [.2 .2 1.5 1.5], 'outerposition', [0 0 1 1]);
cb= colorbar; set(cb, 'location', 'southoutside', 'fontsize', 22, 'position', [.73 .12 .26 .06]); title(cb, ' (uVm ^2)');

% ylims=[-5 5];
% f1=figure;
% % title(['Topo RT Spline ' num2str(spline) ' HPF ' labels{HPF} '_BL' num2str(BL)], 'fontsize', 22, 'fontweight', 'bold');
% %ylims=[min(nanmean(nanmean(TopoEv,2),3))-1 max(nanmean(nanmean(TopoEv,2),3))+1]; ylims= [-20 20];
% h=topoplot(nanmean(nanmean(TopoDiff(:,4,subs),2),3), chanlocs,'plotchans',[1:128],'maplimits', ylims, ...
% 'electrodes','off', 'whitebk', 'off', 'headrad',.65, 'plotrad', .71, 'emarker2', {[3 4],'o','k',10,1}, 'hcolor', 'k', 'conv', 'on');
% %set(h, 'position', [.2 .2 1.5 1.5], 'outerposition', [0 0 1 1]);
% cb= colorbar; set(cb, 'location', 'southoutside', 'fontsize', 22, 'position', [.73 .12 .26 .06]); title(cb, ' (uVm ^2)');

%set(f1,'position',scrsz);
%set(f1,'position',scrsz);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 10])

clear CPPchans; 
CPPChans = [1 2 3 4 5 18 19 20 21 31 32];
for s=1:length(allsubj)
    for d=1:5
        if s~=19
            if s==11 && d==5
            dd=4;
        else dd=d;
            end
        
    thisdat=nanmean(TopoRT(CPPChans,dd,s),2);
%     zdat=zscore(thisdat+10);
%     outliers=zdat<-3 | zdat>3; 
%     thisdat(outliers)=nan; 
    [m idx]=sort(thisdat,'descend');
    choose=sum(isnan(m))+1; % avoid choosing a channel with nans 
    CPPchans{s}(d,:)= CPPChans(idx(choose:choose+3));
        end
    end
end


TopoData_rlock=squeeze(nanmean(PlotDataTraining.ERPr,4));
TopoData_slock=squeeze(nanmean(PlotDataTraining.ERP,4));


%% Load Data
Delays=[800 1200 1600];

load([folder 'LPF_HPF_0_05\SingleTrial_CPP__HPF_0_05_CSD4']);

Delay(:,2)=Delay(:,2)*(1000/512);
RT2=RT+Delay(:,2)';

%% Exclude Artifact/Blink Trials 
AnyBlinks = (blinks.Epoch)>0;
AnyArts = (artifs.Epoch)>0;
Exclude = AnyBlinks>0 | sum(AnyArts,1)>0;

CPP(:,:,Exclude)=nan;  CPPr(:,:,Exclude)=nan;

CPP(:,:,SubID==1 & Sessid==1)=nan; CPPr(:,:,SubID==1 & Sessid==1)=nan; 

%% Choose Baseline Correction Window
BL1=slock_wind{3}>-300 & slock_wind{3}<-200;
BL2=slock_wind{3}>500 & slock_wind{3}<550; 
BL=BL2; 

CPPBL=nanmean(CPP(:,BL,:),2); 
CPP=CPP-repmat(CPPBL, [1 size(CPP,2) 1]); 
CPPr=CPPr-repmat(CPPBL, [1 size(CPPr,2) 1]); 


ZCPP=nan(1,size(CPP,3));ZCPP2=nan(1,size(CPP,3));

%% Condense Matrices 
evwind=[-50 50]; rwin=[-150 -50];
for s=1:length(allsubj)
    if s~=19
    for d=sessions{s}
        
        idx=SubID==s & Sessid==d;
        idxch=ismember(CPPch, CPPchans{s}(d,:));
        CPP2(:,SubID==s & Sessid==d)=nanmean(CPP(idxch,:,idx),1);
        CPPr2(:,SubID==s & Sessid==d)=nanmean(CPPr(idxch,:,idx),1);
        thisdat=[]; thisdat=squeeze(nanmean(CPPr(idxch, rlock_wind>rwin(1) & rlock_wind<rwin(2), idx),1));
        CPPRT(idx)=nanmean(filtfilt(B,A, double(thisdat)),1);
        for z=1:3
            idx2=SubID==s & Sessid==d & Delay(:,1)'==z;
            timeidx=slock_wind{z}>Delays(z)+evwind(1) & slock_wind{z}<Delays(z)+evwind(2);
            thisdat=[]; thisdat=squeeze(nanmean(CPP(idxch,timeidx,idx2),1));
            CPPEvOn(idx2)=nanmean(filtfilt(B,A, double(thisdat)),1);
        end
    end
    end
end

checknans=isnan(CPPRT);
clear CPP CPPr
save(['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\Analyses\Temporal_Uncertainty_Delays\PlottingData\CPP_PlottingData'], 'Exclude', 'CPP*', 'slock_wind', 'rlock_wind',  'Topo*','SubID', 'Acc', 'Block', 'Delay','Resp', 'RT', 'RT2', 'Sessid', 'Targ', 'subs', 'colours*')






