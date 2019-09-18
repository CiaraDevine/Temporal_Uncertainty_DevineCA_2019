%% 
close all
clear all
load('C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\Analyses\Temporal_Uncertainty_Delays\PlottingData\CPP_PlottingData'); 
delays=[800 1200 1600];
chanlocs = readlocs('cap128.loc');
fs=512;
numsubs=22;
subs=[2:5 7:18 20:22];
colours_sessions = {[1, .65, .67] [.95, .35, .40] [1, 0, .06] [.5, 0, 0] [0, 0, 0]}; % Red theme colours


%% Make 6hz filter for visualising CPP
[B,A]=butter(4,(6*2)/fs); % creating filter for ERP waveforms
%% RT Binned  CPP
ZCPP=nan(1, size(CPPRT,2));
CheckNans=isnan(CPPRT);
slidewinds=[800:50:1600];
FPs=[1 2 3];
clear Plot* RTBins*
for s=1:numsubs
    if s~=19
    
        
        for z=1:3
            PlotCPPEv1(z,s)=nanmean(nanmean(CPP2(slock_wind{z}>delays(z)-50 & slock_wind{z}<delays(z)+50, SubID==s & ismember(Delay(:,1)', z))));
            numbins=6;
            RTBins(:,z,s)=prctile(RT2(ismember(Delay(:,1)', z) & SubID==s & ~CheckNans), linspace(0,100,numbins+1));
            for b=1:numbins
                clear idx*;
                idxRTtemp=RT2>RTBins(b,z,s) & RT2<RTBins(b+1,z,s) & ismember(Delay(:,1)', z) & SubID==s & ~CheckNans;
                Zdata=[]; Zdata=CPPRT(idxRTtemp);
                ZCPP(idxRTtemp)= zscore(Zdata);
                
                idxRT=idxRTtemp & ZCPP>-3 & ZCPP<3;
                idxOutliers=idxRTtemp & ZCPP<-3 | ZCPP>3;
                Plotnumtrials(b,z,s)=sum(idxRT);
                PlotnumOutliers(b,z,s)=sum(idxOutliers);
                PlotRT2(b,z,s)=nanmean(RT2(idxRT));
                PlotAcc(b,z,s)=(nansum(idxRT & Acc==1)/nansum(idxRT))*100;
                PlotCPP(:,b,z,s)=squeeze(nanmean(CPP2(:, idxRT),2));
                PlotCPPr(:,b,z,s) =squeeze(nanmean(CPPr2(:, idxRT),2));
                PlotCPPEv(b,z,s)=nanmean(CPPEvOn(idxRTtemp));
                PlotCPPEv_ttest(b,z,s)=nanmean(nanmean(CPP2(slock_wind{z}>delays(z)-50 & slock_wind{z}<delays(z)+50, idxRTtemp)));
                PlotCPPRT(b,z,s) =nanmean(CPPRT(idxRTtemp));

                
                idx=slock_wind{z}>delays(z)-250 & slock_wind{z}<delays(z)+50;
                x=slock_wind{z}(idx)'; y=filtfilt(B,A, double(squeeze(nanmean(CPP2(idx,idxRTtemp),2))));
                P=polyfit(x,y,1);
                PreEvSlope(b,z,s)=P(1);
                idx=slock_wind{z}>delays(z)+200 & slock_wind{z}<delays(z)+600;
                x=slock_wind{z}(idx)'; y=filtfilt(B,A, double(squeeze(nanmean(CPP2(idx,idxRTtemp),2))));
                P=polyfit(x,y,1);
                PostEvSlope(b,z,s)=P(1);
                idx=rlock_wind>-500 & rlock_wind<-200;
                x=rlock_wind(idx)'; y=filtfilt(B,A, double(squeeze(nanmean(CPPr2(idx,idxRTtemp),2))));
                P=polyfit(x,y,1);
                PreRTSlope(b,z,s)=P(1);
            end

        end
    else
    end
end

%% 
 ploton=1; 
for k =ploton
figure; hold on; 

ylims=[-2.99 22.99];

sp(3)= subplot(2,5,3); hold on;
for z=1:3
    Pl(z)=plot(rlock_wind, filtfilt(B,A, double(nanmean(nanmean(PlotCPPr(:,:,z,subs),2),4))) ,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    xlim([-810 400]); ylim(ylims);
    xlabel('Time (ms)', 'fontsize', 18); %ylabel('Power', 'fontsize', 18);
    line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 18, 'ycolor', 'w', 'yticklabels', [], 'ytick', [],'xtick', [-800 -400 0 400], 'xticklabels', {'-800' '-400' '0' '400'});
end


sp(2)=subplot(2,5,2); hold on;
%title(['CPP ' labels{HPF} ' CSD SPline ' num2str(spline)], 'fontsize', 20, 'fontweight', 'bold');
for z=1:3
    Pl(z)=plot(slock_wind{z}, filtfilt(B,A, double(nanmean(nanmean(PlotCPP(1:length(slock_wind{z}),:,z,subs),2),4))) ,'linewidth', 4, 'linestyle', '-',  'color', colours_sessions{z+2});
    %Amp=nanmean(temp(slock_wind{z}>delays(z) & slock_wind{z}<delays(z)+25));
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
end
xlim([500 3600]); ylim(ylims);
line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
set(gca, 'fontsize', 18, 'xtick', [800: 400 : 3600], 'xticklabels', {'800' '1200' '1600' '2000' '2400' '2800' '3200' ''});
%ylabel('CPP Amplitude uV', 'fontsize', 26);
ylabel('CPP [uV/m ^2]', 'fontsize', 24);
xlabel('Time (ms)', 'fontsize', 22);

%lg=legend([Pl(1) Pl(2) Pl(3)], '800ms FP', '1200ms FP', '1600ms FP');
%set(lg, 'fontsize', 24, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast', 'box', 'off');
% h=axes; set(h,'position', [0 .5 .49 .49], 'visible', 'off');label=text(.02, .97, 'A'); set(label,'fontsize', 34, 'fontweight', 'bold');
% h1=axes; set(h1,'position', [.51 .5 .49 .49], 'visible', 'off');label=text(.3, .97, 'B'); set(label,'fontsize', 34, 'fontweight', 'bold');


ylims=[-15 15];
sp(1)=subplot(2,5,1); topoplot(nanmean(nanmean(TopoRT(:,:,subs),2),3), chanlocs,'plotchans',[1:128],'maplimits', ylims, ...
    'electrodes','off', 'whitebk', 'off', 'headrad',.65, 'plotrad', .66, 'emarker2', {[],'o','k',10,1}, 'hcolor', 'k', 'conv', 'on');
%set(h, 'position', [.2 .2 1.5 1.5], 'outerposition', [0 0 1 1]);
cb= colorbar; set(cb, 'location', 'southoutside', 'xtick', [-10 0 10], 'fontsize', 18, 'position', [.015 .48 .1 .01]); title(cb, 'uV/m^2 ', 'fontsize', 18);

ylims=[-2.99 22.99];
sp(4)=subplot(2,5,4); hold on; box off; tt=title ('Pre RT'); 
for z=1:3
    %pl(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(PlotCPPEv(:,z,subs),3), 'linewidth',1.5, 'color', colours_sessions{z+2}); hold on;
    pl2(z)=plot(nanmean(PlotRT2(:,z,subs),3),nanmean(PlotCPPRT(:,z,subs),3),'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',1, 'color', colours_sessions{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [], 'xticklabel', {}, 'yticklabel', {});
ylabel(' ','fontsize',19); %xlabel('RT (ms)', 'fontsize',20);
% lg=legend([pl(3) pl2(3)],'At Response', 'Pre Contrast Change');
% set(lg, 'fontsize', 24, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast');
set(tt, 'fontsize', 22); 


sp(5)=subplot(2,5,5); hold on; box off; tt=title ('Pre Evidence'); 
for z=1:3
    pl(z)=plot(nanmean(PlotRT2(:,z,subs),3),nanmean(PlotCPPEv(:,z,subs),3), 'linewidth',2, 'color', colours_sessions{z+2}); hold on;
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(PlotCPPRT(:,z,subs),3),'linestyle', '-', 'linewidth',3, 'color', colours_sessions{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',1, 'color', colours_sessions{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [], 'xticklabel', {});
ylabel('CPP [uVm^2] ','fontsize',19); %xlabel('RT (ms)', 'fontsize',20);
set(tt, 'fontsize', 22); 


ylims=[-.01 .045];
sp(6)=subplot(2,5,6); hold on; box off; 
for z=1:3
    %pl(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(PlotCPPEv(:,z,subs),3), 'linewidth',1.5, 'color', colours_sessions{z+2}); hold on;
    pl2(z)=plot(nanmean(PlotRT2(:,z,subs),3), nanmean(PreRTSlope(:,z,subs),3),'linestyle', '-', 'linewidth',2, 'color', colours_sessions{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',1, 'color', colours_sessions{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [800:400:3200], 'xticklabel', {'800' '' '1600' '' '2400' '' '3200' },...
    'ytick',[0 .02 .04], 'yticklabel', {});
ylabel('','fontsize',19); 
xlabel('RT (ms)', 'fontsize',20);
% lg=legend([pl(3) pl2(3)],'At Response', 'Pre Contrast Change');
% set(lg, 'fontsize', 24, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast');


sp(7)=subplot(2,5,7); hold on; box off; 
for z=1:3
    pl(z)=plot(nanmean(PlotRT2(:,z,subs),3),nanmean(PreEvSlope(:,z,subs),3), 'linewidth',2, 'color', colours_sessions{z+2}); hold on;
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(PlotCPPRT(:,z,subs),3),'linestyle', '-', 'linewidth',3, 'color', colours_sessions{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',1, 'color', colours_sessions{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [800:400:3200], 'xticklabel', {'800' '' '1600' '' '2400' '' '3200' },...
    'ytick',[0 .02 .04], 'yticklabel', {'0' '.02' '.04'});
ylabel('CPP Slope','fontsize',19); 
xlabel('RT (ms)', 'fontsize',20);




set(sp(1), 'position', [-.059 .49 .25 .25]);
set(sp(2), 'position', [.21 .4 .3 .58]);
set(sp(3), 'position', [.53 .4 .11 .58]);
set(sp(5), 'position', [.735 .70 .119 .21]);
set(sp(4), 'position', [.88 .70 .119 .21]);
set(sp(7), 'position', [.735 .43 .118 .21]);
set(sp(6), 'position', [.88 .43 .118 .21]);


h=axes; set(h,'position', [0 .01 .99 .99], 'visible', 'off');% Create a new axis that is invisible. Otherwise all text will be positioned relevant to your most recent lot which you might not want.
text('String', 'Decision Formation','fontsize', 30, 'fontweight', 'bold', 'position',[.001, .94]);
text('String', 'A','fontsize', 34, 'fontweight', 'bold', 'position',[.001, .79]);
text('String', 'B','fontsize', 34, 'fontweight', 'bold', 'position',[.17, .97]);
text('String', 'C', 'fontsize', 34, 'fontweight', 'bold','position',[.515, .97]);
text('String', 'D','fontsize', 28, 'fontweight', 'bold','position',[.68, .92]);
text('String','E','fontsize', 28, 'fontweight', 'bold','position',[.865, .92]);
text('String', 'F','fontsize', 28, 'fontweight', 'bold','position',[.68, .66]);
text('String','G','fontsize', 28, 'fontweight', 'bold','position',[.865, .66]);
end %%



%%



figure; hold on; 
set(gca, 'fontsize', 22);
for z=1:3
    pl(z)=plot(nanmean(PlotRT2(:,z,subs),3), nanmean(PreRTSlope(:,z,subs),3), 'linestyle', '-','linewidth',4, 'color', colours_sessions{z+2});
    pl(z)=plot(nanmean(PlotRT2(:,z,subs),3), nanmean(PreEvSlope(:,z,subs),3), 'linestyle', '-','linewidth',2, 'color', colours_sessions{z+2});
   % pl2(d)=plot(nanmean(PlotRT2(:,z,subs),3), nanmean(PlotCPPEv(:,z,subs),3), 'linestyle', '-','linewidth',2, 'color', colours_sessions{z+2});
    set(gca, 'fontsize', 22, 'xtick', [1600 2000 2400 2800]); 
end
xlim([1550 2950]); ylim([0 .05]);
ylabel('CPP uV/m^2', 'fontsize', 24); 
xlabel('RT(ms)', 'fontsize', 24); 
% lg=legend([pl(1) pl(2) pl(3) pl(4) pl(5) pl1(5) pl2(5)], 'Day 1','Day 2','Day 3','Day 4','Day 5', '@ RT', '@ Evidence Onset'); 
% set(lg,'fontsize', 26, 'box', 'off'); 


Inputdata=PlotCPPEv(:,:,subs);
ANOVA.CPPAmpEv=RMAOV(Inputdata);
ANOVAContrasts.CPPAmpEv=RMAOV(nanmean(Inputdata), [],[],[], 1);
ANOVAContrasts.CPPAmpEv_RT=RMAOV(nanmean(Inputdata(:,3,:),2), [],[],[], 1);
figure; hist(ANOVA.CPPAmpEv.Normality.StdResiduals2);

Inputdata=PlotCPPRT(:,:,subs);
ANOVA.CPPAmpRT=RMAOV(Inputdata);
ANOVAContrasts.CPPAmpRT=RMAOV(nanmean(Inputdata), [],[],[], 1);
%figure; hist(ANOVA.CPPAmpRT.Normality.StdResiduals2);

Inputdata=PreEvSlope(:,:,subs);
ANOVA.CPPSlopePreEv=RMAOV(Inputdata);
ANOVAContrasts.CPPSlopePreEv=RMAOV(nanmean(Inputdata), [],[],[], 1);
[hsl psl csl tsl]=ttest(nanmean(PreEvSlope(:,2,subs)), nanmean(PreEvSlope(:,3,subs)));

Inputdata=PreRTSlope(:,:,subs);
ANOVA.CPPSlopePreRT=RMAOV(Inputdata);
ANOVAContrasts.CPPSlopePreRT=RMAOV(nanmean(Inputdata), [],[],[], 1);
%figure; hist(ANOVA.CPPSlopePreRT.Normality.StdResiduals2);

% Test CPP Amplitude at 800, 1200 and 1600ms against 0
[h p c t]=ttest(nanmean(PlotCPPEv1(1,subs),1),0);
[h2 p2 c2 t2]=ttest(squeeze(nanmean(PlotCPPEv1(2,subs),1)),0);
[h3 p3 c3 t3]=ttest(squeeze(nanmean(PlotCPPEv1(3,subs),1)),0);
figure; hist(squeeze(nanmean(PlotCPPEv1(3,subs),1)));

% Test CPP Slope prior to 800, 1200 and 1600ms against 0
[h pp cc tt]=ttest(nanmean(PreEvSlope(:,1,subs)),0);
[h2 pp2 cc2 tt2]=ttest(nanmean(PreEvSlope(:,2,subs)),0);
[h3 pp3 cc3 tt3]=ttest(nanmean(PreEvSlope(:,3,subs)),0);
figure; hist(squeeze(nanmean(PreEvSlope(:,3,subs),1)));



