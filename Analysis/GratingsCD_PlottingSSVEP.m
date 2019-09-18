%% Plotting Mu/Beta - Delays 
close all; 
clear all;
load('C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\Analyses\Temporal_Uncertainty_Delays\PlottingData\SSVEP_PlottingData'); 
delays=[800 1200 1600];
chanlocs = readlocs('cap128.loc');
colours_red = {[1, .65, .67] [.95, .35, .40] [1, 0, .06] [.5, 0, 0] [0, 0, 0]}; % Red theme colours

numsub=22;
subs=[2:5 7:18 20:22];

%% 
SSVEPDiff(:, Resp==12)=SSVEPSNR_20Hz(:,Resp==12)-SSVEPSNR_25Hz(:,Resp==12);
SSVEPDiff(:, Resp==13)=SSVEPSNR_25Hz(:,Resp==13)-SSVEPSNR_25Hz(:,Resp==13);
SSVEPDiff_Cum=cumsum(SSVEPDiff,1);
%%
windR=TimePointsR>-300 & TimePointsR<-200; 
preevslopewin=[-300 0];
preRTslopewin=[-500 -200];


ZSSVEP=nan(1, size(SubID,2));
for s=1:numsub
    for z=1:3
       idxLt= Targ==8 & Delay(1,:)==z & SubID==s; idxRt= Targ==9 & Delay(1,:)==z & SubID==s & RT>500; 
        idxLr= Resp==12 & Delay(1,:)==z & SubID==s; idxRr= Resp==13 & Delay(1,:)==z & SubID==s & RT>500; 
        ZSSVEP(idxLt)=zscore(nanmean(SSVEPSNR_20Hz(windR, idxLt)));
        ZSSVEP(idxRt)=zscore(nanmean(SSVEPSNR_25Hz(windR, idxRt)));
        Outliers=ZSSVEP>3 | ZSSVEP<-3;
        
        SSVEP20Lt(:,z,s)=nanmean(SSVEPSNR_20Hz(:,idxLt & ~Outliers),2);
        SSVEP20Rt(:,z,s)=nanmean(SSVEPSNR_20Hz(:,idxRt & ~Outliers),2);
        SSVEP25Lt(:,z,s)=nanmean(SSVEPSNR_25Hz(:,idxLt & ~Outliers),2);
        SSVEP25Rt(:,z,s)=nanmean(SSVEPSNR_25Hz(:,idxRt & ~Outliers),2);
        SSVEPr20Lt(:,z,s)=nanmean(SSVEPSNRr_20Hz(:,idxLt & ~Outliers),2);
        SSVEPr20Rt(:,z,s)=nanmean(SSVEPSNRr_20Hz(:,idxRt & ~Outliers),2);
        SSVEPr25Lt(:,z,s)=nanmean(SSVEPSNRr_25Hz(:,idxLt & ~Outliers),2);
        SSVEPr25Rt(:,z,s)=nanmean(SSVEPSNRr_25Hz(:,idxRt & ~Outliers),2);
        if z==3
        SSVEPDiff_Fast(:,s)=nanmean(SSVEPDiff(:,RT<150 & Delay(1,:)==z & SubID==s),2);    
        SSVEPCumDiff_Fast(:,s)=nanmean(SSVEPDiff_Cum(:,RT<150 & Delay(1,:)==z & SubID==s),2);    
        FastNumTrials(s)=sum(RT<250 & Delay(1,:)==z & SubID==s);
        end
        numbins=8;
        RTBinsLt(:,z,s)=prctile(RT2(idxLt & RT<1800), linspace(0,100,numbins+1)); RTBinsRt(:,z,s)=prctile(RT2(idxRt & RT<1800), linspace(0,100,numbins+1));
        RTBinsLr(:,z,s)=prctile(RT2(idxLr & RT<1800), linspace(0,100,numbins+1)); RTBinsRr(:,z,s)=prctile(RT2(idxRr & RT<1800), linspace(0,100,numbins+1));
        for b=1:numbins
            
            idxRTLt = RT2>RTBinsLt(b,z,s) & RT2<RTBinsLt(b+1,z,s) & idxLt; numtrialsLt(b,z,s)=sum(idxRTLt);
            idxRTRt = RT2>RTBinsRt(b,z,s) & RT2<RTBinsRt(b+1,z,s) & idxRt; numtrialsRt(b,z,s)=sum(idxRTRt);
            idxRTLr = RT2>RTBinsLr(b,z,s) & RT2<RTBinsLr(b+1,z,s) & idxLr; numtrialsLr(b,z,s)=sum(idxRTLr);
            idxRTRr = RT2>RTBinsRr(b,z,s) & RT2<RTBinsRr(b+1,z,s) & idxRr; numtrialsRr(b,z,s)=sum(idxRTRr);
            RTBins(b,z,s)=nanmean(RT2(idxRTLt | idxRTRt));
           
            
            SSVEP20Lt2(:,b,z,s) = nanmean(SSVEPSNR_20Hz(:,idxRTLt & ~Outliers),2);
            SSVEP20Rt2(:,b,z,s) = nanmean(SSVEPSNR_20Hz(:,idxRTRt & ~Outliers),2);
            SSVEP25Lt2(:,b,z,s) = nanmean(SSVEPSNR_25Hz(:,idxRTLt & ~Outliers),2);
            SSVEP25Rt2(:,b,z,s) = nanmean(SSVEPSNR_25Hz(:,idxRTRt & ~Outliers),2);
            
            SSVEPr20Lt2(:,b,z,s) = nanmean(SSVEPSNRr_20Hz(:,idxRTLt & ~Outliers),2);
            SSVEPr20Rt2(:,b,z,s) = nanmean(SSVEPSNRr_20Hz(:,idxRTRt & ~Outliers),2);
            SSVEPr25Lt2(:,b,z,s) = nanmean(SSVEPSNRr_25Hz(:,idxRTLt & ~Outliers),2);
            SSVEPr25Rt2(:,b,z,s) = nanmean(SSVEPSNRr_25Hz(:,idxRTRt & ~Outliers),2);
            
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


figure; hold on; 
xidx=TimePoints<(delays(3)+1800);
plot(TimePoints(xidx), nanmean(SSVEPDiff(xidx,subs),2), 'linewidth', 3, 'color', 'k'); 
xlim([0 3400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
line([delays(3) delays(3)],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
set(gca, 'fontsize', 20);

figure; hold on; 
xidx=TimePoints<(delays(3)+1800);
plot(TimePoints(xidx), nanmean(SSVEPCumDiff_Fast(xidx,subs),2), 'linewidth', 3, 'color', 'k'); 
xlim([0 3400]); %ylim(ylims); %ylabel('Power', 'fontsize', 18);
line([delays(3) delays(3)],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
set(gca, 'fontsize', 20);
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
