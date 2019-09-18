%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

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

subs=[2:5 7:18 20:22];
%% Make 6hz filter for visualising N2_
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
for z=1:3 
    TopoEv(:,z,d,:) = nanmean(nanmean(PlotDataTraining.ERP(:,slock_wind{1}>delays(z)+200 & slock_wind{1}<delays(z)+300,d, :, z,:),2),4);
end
end

%ylims=[min(nanmean(nanmean(TopoRT,2),3))-3 max(nanmean(nanmean(TopoRT,2),3))+3]; %

ylims=[-10 10];
f1=figure;
% title(['Topo RT Spline ' num2str(spline) ' HPF ' labels{HPF} '_BL' num2str(BL)], 'fontsize', 22, 'fontweight', 'bold');
%ylims=[min(nanmean(nanmean(TopoEv,2),3))-1 max(nanmean(nanmean(TopoEv,2),3))+1]; ylims= [-20 20];
h=topoplot(nanmean(nanmean(TopoEv(:,1,subs),2),3), chanlocs,'plotchans',[1:128],'maplimits', ylims, ...
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

clear N2_chans; 
N2_Chans1 = [6:9 125:126];
N2_Chans2 = [ 35:38 44 45];
for s=1:length(allsubj)
    for d=1:5
        if s~=19
            if s==11 && d==5
            dd=4;
        else dd=d;
            end
        
    thisdat1=nanmean(TopoEv(N2_Chans1,:,dd,s),2);
    thisdat2=nanmean(TopoEv(N2_Chans2,:,dd,s),2);
%     zdat=zscore(thisdat+10);
%     outliers=zdat<-3 | zdat>3; 
%     thisdat(outliers)=nan; 
    [m1 idx1]=sort(thisdat1,'ascend');
    [m2 idx2]=sort(thisdat2,'ascend');
    choose1=sum(isnan(m1))+1; % avoid choosing a channel with nans 
    choose2=sum(isnan(m2))+1; % avoid choosing a channel with nans 
    N2_chans{s}(d,1,:)= N2_Chans1(idx1(choose1:choose1+2));
    N2_chans{s}(d,2,:)= N2_Chans2(idx2(choose2:choose2+2));
        end
    end
end


%% Load Data
Delays=[800 1200 1600];

clear PlotData*
load([folder 'LPF_HPF_0_05\SingleTrial_N2___HPF_0_05_CSD4']);
%load(['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\GratingsDelaysData\ERPMatrices\LPF_HPF_0_05\SingleTrial_CPP__HPF_0_05_CSD4']);
%load([plotfolder, 'N2_ChansSessions']); 

Delay(:,2)=Delay(:,2)*(1000/512);
RT2=RT+Delay(:,2)';

%% Exclude Artifact/Blink Trials 
AnyBlinks = (blinks.Epoch)>0;
AnyArts = (artifs.Epoch)>0;
Exclude = AnyBlinks>0 | sum(AnyArts,1)>0;

N2(:,:,Exclude)=nan;  N2r(:,:,Exclude)=nan;

N2(:,:,SubID==1 & Sessid==1)=nan; N2r(:,:,SubID==1 & Sessid==1)=nan; 

N2_BL=nanmean(N2(:,1:50,:),2); 
N2=N2-repmat(N2_BL, [1 size(N2,2) 1]); 
N2r=N2r-repmat(N2_BL, [1 size(N2r,2) 1]); 


%ZN2_=nan(1,size(N2_,3)); ZN2_2=nan(1,size(N2_,3));

%% Condense Matrices 
evwind=[-50 50]; rwin=[-150 -50];
for s=1:length(allsubj)
    if s~=19
        for d=sessions{s}
            for z=1:3
                idx=SubID==s & Sessid==d & Delay(:,1)'==z;
                idxch=find([ismember(N2_ch, N2_chans{s}(d,1,:)) + ismember(N2_ch, N2_chans{s}(d,2,:))]);
                N2_2(:,idx)=nanmean(N2(idxch,:,idx),1);
            end
        end
    end
end
wind=(-50:400)*1.953;
N2(:,:,Exclude)=nan;  N2r(:,:,Exclude)=nan;
clear N2 N2r
save(['C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\Analyses\Temporal_Uncertainty_Delays\PlottingData\N2_PlottingData'], 'Exclude', 'N2*', 'slock_wind', 'rlock_wind',  'wind','Topo*','SubID', 'Acc', 'Block', 'Delay','Resp', 'RT', 'RT2', 'Sessid', 'Targ', 'subs', 'colours*')



%% Waveforms N2_
% RT Binned  N2_
clear Plot* RTBins*
for s=1:length(allsubj)
    if s~=19
        for z=1:3
            PlotN2(:,z,s)=squeeze(nanmean(N2_2(:, SubID==s & ismember(Delay(:,1)', z)),2));
            N2Amp(z,s)=nanmean(nanmean(N2_2(128:150, SubID==s & ismember(Delay(:,1)', z)),2),1);
        end
    else
    end
end

%% 
 subs=[2:5 7:18 20:22];
 
%% N2 Stim-locked 
figure; hold on; 

ylims2=[-7.99 7.99];
sp(1)=subplot(1,2,1); topoplot(nanmean(nanmean(TopoEv(:,:,subs),2),3), chanlocs,'plotchans',[1:128],'maplimits', ylims2, ...
    'electrodes','off', 'whitebk', 'off', 'headrad',.65, 'plotrad', .66, 'emarker2', {[],'o','k',10,1}, 'hcolor', 'k', 'conv', 'on');
%set(h, 'position', [.2 .2 1.5 1.5], 'outerposition', [0 0 1 1]);
%cb= colorbar; set(cb, 'location', 'eastoutside', 'fontsize', 22, 'position', [.95 .63 .015 .25]); %title(cb, ' (uVm ^2)');

ylims=[-12.99 4.99]; 
sp(2)=subplot(1,2,2); hold on;
%title(['N2_ ' labels{HPF} ' CSD SPline ' num2str(spline)], 'fontsize', 20, 'fontweight', 'bold');
for z=1:3
    wind=(-50:400)*1.953;
    Pl(z)=plot(wind, filtfilt(B,A, double(nanmean(nanmean(PlotN2(1:length(wind),z,subs),2),3))) ,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    %Amp=nanmean(temp(slock_wind{z}>delays(z) & slock_wind{z}<delays(z)+25));
    %l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    xlabel('Time (ms)', 'fontsize', 22);
end
xlim([-110 810]); ylim(ylims);
line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
set(gca, 'fontsize', 22, 'xtick', [-100:50:800], 'xticklabels', {'-100' '' '' '' '100' '' '200' '' '300' '' '400' '' '500' '' '600' '' '700' '' '800'});
%ylabel('N2_ Amplitude uV', 'fontsize', 26);
ylabel('N2 [uV/m ^2]', 'fontsize', 24);


set(sp(1), 'position', [-.04 .5 .25 .25]);
set(sp(2), 'position', [.26 .4 .4 .58]);

h=axes; set(h,'position', [0 .01 .99 .99], 'visible', 'off');% Create a new axis that is invisible. Otherwise all text will be positioned relevant to your most recent lot which you might not want.
text('String', 'Target Selection','fontsize', 30, 'fontweight', 'bold', 'position',[.01, .94]);
text('String', 'A','fontsize', 34, 'fontweight', 'bold', 'position',[.01, .78]);
text('String', 'B','fontsize', 34, 'fontweight', 'bold', 'position',[.2, .97]);




%% Stats 
[h1 p1 t1 CI1]=ttest(N2Amp(1,subs), 0); 
[h2 p2 t2 CI2]=ttest(N2Amp(2,subs), 0); 
[h3 p3 t3 CI3]=ttest(N2Amp(3,subs), 0); 

ANOVA.N2=RMAOV(N2Amp(:,subs));
