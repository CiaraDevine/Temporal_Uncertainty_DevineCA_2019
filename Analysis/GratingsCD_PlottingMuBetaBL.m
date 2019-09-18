%% Plotting Mu/Beta - Delays 
close all; 
clear all;
load('C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\Analyses\Temporal_Uncertainty_Delays\PlottingData\BetaBL_PlottingData'); 
delays=[800 1200 1600];
chanlocs = readlocs('cap128.loc');
colours_red = {[1, .65, .67] [.95, .35, .40] [1, 0, .06] [.5, 0, 0] [0, 0, 0]}; % Red theme colours

numsub=22;
subs=[2:5 7:18 20:22];

%%

preevslopewin=[-300 -0];
preRTslopewin=[-500 -200];
checknans=isnan(BetaRT(1,:));
ZFFT=nan(1, size(SubID,2));
for s=1:numsub
    for z=1:3
        idx=Delay(1,:)==z &  SubID==s;% & RT>1000 & RT<1500; 
        
        BetacontraLHtemp = nanmean(Beta(1,:, Resp==13 & idx),3);
        BetacontraRHtemp = nanmean(Beta(2,:, Resp==12 & idx),3);
        BetaContra(:,z,s)= (BetacontraLHtemp + BetacontraRHtemp)/2; 
        BetaipsiLHtemp = nanmean(Beta(1,:, Resp==12 & idx),3);
        BetaipsiRHtemp = nanmean(Beta(2, :, Resp==13 & idx),3);
        BetaIpsi(:,z,s)= (BetaipsiLHtemp + BetaipsiRHtemp)/2;
        
        BetaRcontraLHtemp = nanmean(BetaR(1,:, Resp==13 & idx),3);
        BetaRcontraRHtemp = nanmean(BetaR(2,:, Resp==12 & idx),3);
        BetaRContra(:,z,s)= (BetaRcontraLHtemp + BetaRcontraRHtemp)/2; 
        BetaRipsiLHtemp = nanmean(BetaR(1,:, Resp==12 & idx),3);
        BetaRipsiRHtemp = nanmean(BetaR(2,:,Resp==13 & idx),3);
        BetaRIpsi(:,z,s)= (BetaRipsiLHtemp + BetaRipsiRHtemp)/2;
        
        
        numbins=6;
        RTBins(:,z,s)=prctile(RT2(Delay(1,:)==z & SubID==s & RT<1800), linspace(0,100,numbins+1));
        for b=1:numbins
            idxRTtemp = RT2>RTBins(b,z,s) & RT2<RTBins(b+1,z,s) & idx & ~checknans;
            Plotnumtrials(b,z,s)=sum(idxRTtemp);
            Plot_RTBins(b,z,s)=nanmean(RT2(idxRTtemp));
            ZFFT(idxRTtemp)=squeeze(zscore(nanmean(BetaRT(:, idxRTtemp),1)));
            idxRT = RT2>RTBins(b,z,s) & RT2<RTBins(b+1,z,s) & idx & ~checknans;% & ZFFT>-3 & ZFFT<3;
            
            BetaContraLH(b,z,s) = nanmean(BetaEv(1, Resp==13 & idxRT),2); 
            BetaContraRH(b,z,s) = nanmean(BetaEv(2, Resp==12 & idxRT),2); 
            BetaContraBins(b,z,s) = (BetaContraLH(b,z,s) + BetaContraRH(b,z,s))/2;
 

            BetaIpsiLH(b,z,s) = nanmean(BetaEv(1, Resp==12 & idxRT),2);
            BetaIpsiRH(b,z,s) = nanmean(BetaEv(2, Resp==13 & idxRT),2); 
            BetaIpsiBins(b,z,s) = (BetaIpsiLH(b,z,s) + BetaIpsiRH(b,z,s))/2;
            
            BetaRContraLH(b,z,s) = nanmean(BetaRT(1, Resp==13 & idxRT),2); 
            BetaRContraRH(b,z,s) = nanmean(BetaRT(2, Resp==12 & idxRT),2);
            BetaRContraBins(b,z,s) = (BetaRContraLH(b,z,s) + BetaRContraRH(b,z,s))/2;
           
            BetaRIpsiLH(b,z,s) = nanmean(BetaRT(1, Resp==12 & idxRT),2);
            BetaRIpsiRH(b,z,s) = nanmean(BetaRT(2, Resp==13 & idxRT),2);
            BetaRIpsiBins(b,z,s) = (BetaRIpsiLH(b,z,s) + BetaRIpsiRH(b,z,s))/2;
            
            % Slope Measurements
            BetaContraLHtemp = nanmean(Beta(1,:, Resp==13 & idxRT),3); BetaContraRHtemp = nanmean(Beta(2,:, Resp==12 & idxRT),3);
            BetaContra2(:,b,z,s)= (BetaContraLHtemp + BetaContraRHtemp)/2; 
            BetaIpsiLHtemp = nanmean(Beta(1,:, Resp==12 & idxRT),3); BetaIpsiRHtemp = nanmean(Beta(2, :, Resp==13 & idxRT),3);
            BetaIpsi2(:,b,z,s)= (BetaIpsiLHtemp + BetaIpsiRHtemp)/2;
            BetaRContraLHtemp = nanmean(BetaR(1,:, Resp==13 & idxRT),3); BetaRContraRHtemp = nanmean(BetaR(2,:, Resp==12 & idxRT),3);
            BetaRContra2(:,b,z,s)= (BetaRContraLHtemp + BetaRContraRHtemp)/2; 
            BetaRIpsiLHtemp = nanmean(BetaR(1,:, Resp==12 & idxRT),3); BetaRIpsiRHtemp = nanmean(BetaR(2, :, Resp==13 & idxRT),3);
            BetaRIpsi2(:,b,z,s)= (BetaRIpsiLHtemp + BetaRIpsiRHtemp)/2;
                idxt=TimePoints>delays(z)+ preevslopewin(1) & TimePoints<delays(z)+ preevslopewin(2);
                x=TimePoints(idxt)'; 
                y=squeeze(BetaContra2(idxt,b,z,s)); y2=squeeze(BetaIpsi2(idxt,b,z,s));
                P=polyfit(x,y,1);P2=polyfit(x,y2,1);
                PreEvSlopeContra(b,z,s)=P(1); PreEvSlopeIpsi(b,z,s)=P2(1);
                idxt=TimePointsR>preRTslopewin(1) & TimePointsR<preRTslopewin(2);
                x=TimePointsR(idxt)'; 
                y=squeeze(BetaRContra2(idxt,b,z,s)); y2=squeeze(BetaRIpsi2(idxt,b,z,s));
                P=polyfit(x,y,1);P2=polyfit(x,y2,1);
                PreRTSlopeContra(b,z,s)=P(1);PreRTSlopeIpsi(b,z,s)=P2(1);
        end
    end
end






%% Plot Beta 
ploton=1; 
for h=ploton

figure; 
ylims=[-1.49 .99];
ylims=[-.09 1.73];
sp(8)= subplot(2,7,8); hold on;
for z=1:3
    Plc(z)=plot(TimePointsR(1:end-1), nanmean(BetaRContra(:,z,subs),3)*-1,'linewidth', 4, 'linestyle', '-',  'color', colours_red{z+2});
    Pli(z)=plot(TimePointsR(1:end-1), nanmean(BetaRIpsi(:,z,subs),3)*-1,'linewidth', 4, 'linestyle', '--',  'color', colours_red{z+2});
    xlim([-810 400]); ylim(ylims); %ylabel('Power', 'fontsize', 18);
    line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
    %text([-150], ylims(2)+1.5, 'Response', 'fontsize', 20);
    set(gca, 'fontsize', 20, 'ycolor', 'w', 'yticklabels', [], 'ytick', [],'xtick', [-800 -400 0 400], 'xticklabels', {'-800' '-400' '0' '400'});
end


sp(7)=subplot(2,7,7); hold on;
%title(['CPP ' labels{HPF} ' CSD SPline ' num2str(spline)], 'fontsize', 20, 'fontweight', 'bold');
for z=1:3
    Plc(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), nanmean(BetaContra(TimePoints<(delays(z)+1800),z,subs),3)*-1,'linewidth', 4, 'linestyle', '-',  'color', colours_red{z+2});
    Plc1(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), nanmean(BetaContra(TimePoints<(delays(z)+1800),z,subs),3)*-1,'linewidth', 4, 'linestyle', '-',  'color', colours_red{z+2});
    Pli(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), nanmean(BetaIpsi(TimePoints<(delays(z)+1800),z,subs),3)*-1,'linewidth', 4, 'linestyle', '--',  'color', colours_red{z+2});
    %Amp=nanmean(temp(slock_wind{z}>delays(z) & slock_wind{z}<delays(z)+25));
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_red{z+2});
    
end
xlabel('Time (ms)', 'fontsize', 18);
xlim([500 3600]); ylim(ylims);
line([0 0],ylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
set(gca, 'fontsize', 20, 'xtick', [800: 400 : 3600], 'xticklabels', {'800' '1200' '1600' '2000' '2400' '2800' '3200' ''});
%ylabel('CPP Amplitude uV', 'fontsize', 26);
ylabel('Beta (uV)', 'fontsize', 26);
lg=legend([Plc1(3) Pli(3)], 'Contra', 'Ipsi');
set(lg, 'fontsize', 26, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast', 'box', 'off');
% h=axes; set(h,'position', [0 .5 .49 .49], 'visible', 'off');label=text(.02, .97, 'A'); set(label,'fontsize', 34, 'fontweight', 'bold');
% h1=axes; set(h1,'position', [.51 .5 .49 .49], 'visible', 'off');label=text(.3, .97, 'B'); set(label,'fontsize', 34, 'fontweight', 'bold');

tempylims=[-.25 .6];
sp(13)=subplot(2,7,13); hold on; 
for z=1:3
    Plc(z)=plot(TimePoints(TimePoints<(delays(z)+1800)), (nanmean(BetaContra(TimePoints<(delays(z)+1800),z,subs)-BetaIpsi(TimePoints<(delays(z)+1800),z,subs),3))*-1 ,'linewidth', 2, 'linestyle', '-',  'color', colours_red{z+2});
    %Amp=nanmean(temp(slock_wind{z}>delays(z) & slock_wind{z}<delays(z)+25));
    l(z)=line([delays(z) delays(z)],[tempylims(1) tempylims(2)], 'linestyle', '-', 'linewidth',1, 'color', colours_red{z+2});  
end
xlabel('Time (ms)', 'fontsize', 14);
xlim([500 3600]); ylim(tempylims);
line([0 0],tempylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
%ylabel('CPP Amplitude uV', 'fontsize', 26);
ylabel('Lat Beta (uV)', 'fontsize', 14);
set(gca, 'fontsize', 10, 'xtick', [800: 400 : 3600], 'xticklabels', {'800' '' '1600' '' '2400' '' '3200' ''});


sp(14)=subplot(2,7,14); hold on; 
for z=1:3
    Plc(z)=plot(TimePointsR(1:end-1), (nanmean(BetaRContra(:,z,subs)-BetaRIpsi(:,z,subs),3))*-1,'linewidth', 2, 'linestyle', '-',  'color', colours_red{z+2});
    %Amp=nanmean(temp(slock_wind{z}>delays(z) & slock_wind{z}<delays(z)+25));
    l(z)=line([delays(z) delays(z)],[tempylims(1) tempylims(2)], 'linestyle', '-', 'linewidth',1, 'color', colours_red{z+2});  
end
xlabel('Time (ms)', 'fontsize', 14);
xlim([-810 400]); ylim(tempylims);
line([0 0],tempylims, 'linestyle', '-', 'linewidth',2, 'color', 'k');
%ylabel('CPP Amplitude uV', 'fontsize', 26);
ylabel('', 'fontsize', 14);
set(gca, 'fontsize', 10, 'ytick', [], 'xtick', [-800 -400 0 400], 'xticklabels', {'-800' '-400' '0' '400'});




ylims=[-.5 .5];
sp(6)=subplot(2,7,6); topoplot(nanmean(nanmean(nanmean(nanmean(TopoRT(:,1,:,:,subs)-TopoRT(:,2,:,:,subs),2),3),4),5)*-1, chanlocs,'plotchans',[1:128],'maplimits', ylims, ...
    'electrodes','off', 'whitebk', 'off', 'headrad',.65, 'plotrad', .66, 'emarker2', {[],'o','k',10,1}, 'hcolor', 'k', 'conv', 'on');
%set(h, 'position', [.2 .2 1.5 1.5], 'outerposition', [0 0 1 1]);
%cb= colorbar; set(cb, 'location', 'eastoutside', 'fontsize', 22, 'position', [.95 .63 .015 .25]); %title(cb, ' (uVm ^2)');

ylims=[-.2 1.1];
ylims=[-.09 1.41];
sp(9)=subplot(2,7,9); hold on; box off; tt=title('Pre RT'); 
for z=1:3
    pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),abs(nanmean(BetaRContraBins(:,z,subs),3)-nanmean(BetaRIpsiBins(:,z,subs),3)),'linestyle', '-', 'linewidth',3, 'color', colours_red{z+2});
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(BetaRContraBins(:,z,subs),3)*-1,'linestyle', '-', 'linewidth',3, 'color', colours_red{z+2});
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(BetaRIpsiBins(:,z,subs),3)*-1,'linestyle', '--', 'linewidth',3, 'color', colours_red{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_red{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [500 1000 1500 2000 2500 3000 3500], 'xticklabel', {}, 'ytick', [0 .4 .8 1.2], 'yticklabels', {});
ylabel('','fontsize',24); %xlabel('RT (ms)', 'fontsize',22);
set(tt, 'fontsize', 24);
% lg=legend([pl(3) pl2(3)],'At Response', 'Pre Contrast Change');
% set(lg, 'fontsize', 24, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast');

%ylims=[-.2 1.1];
ylims=[5.41 7.34]; 
ylims=[-.09 1.41];
sp(10)=subplot(2,7,10); hold on; box off; tt=title('Pre Evidence'); 
for z=1:3
    pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),abs(nanmean(BetaContraBins(:,z,subs),3)-nanmean(BetaIpsiBins(:,z,subs),3)),'linestyle', '-', 'linewidth',3, 'color', colours_red{z+2});
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(BetaContraBins(:,z,subs),3)*-1,'linestyle', '-', 'linewidth',3, 'color', colours_red{z+2});
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(BetaIpsiBins(:,z,subs),3)*-1,'linestyle', '--', 'linewidth',3, 'color', colours_red{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_red{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [500 1000 1500 2000 2500 3000 3500], 'xticklabel', { }, 'ytick', [0 .4 .8 1.2]);
ylabel('uV','fontsize',24); %xlabel('RT (ms)', 'fontsize',24);
set(tt, 'fontsize', 24);


ylims=[-.0002 .0026];
%ylims=[4.85 7.15]; 
sp(11)=subplot(2,7,11); hold on; box off; 
for z=1:3
    pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),abs(nanmean(PreRTSlopeContra(:,z,subs),3)-nanmean(PreRTSlopeIpsi(:,z,subs),3)),'linestyle', '-', 'linewidth',3, 'color', colours_red{z+2});
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3), nanmean(PreRTSlopeContra(:,z,subs),3)*-1,'linestyle', '-', 'linewidth',3, 'color', colours_red{z+2});
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3), nanmean(PreRTSlopeIpsi(:,z,subs),3)*-1,'linestyle', '--', 'linewidth',3, 'color', colours_red{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_red{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [500 1000 1500 2000 2500 3000 3500], 'xticklabel', {'' '1000' '' '2000' '' '3000' '' }, ...
    'ytick', [ -.002  -.001 0 .001 .002 .003], 'yticklabels', {});
ylabel('','fontsize',24); xlabel('RT (ms)', 'fontsize',24);
% lg=legend([pl(3) pl2(3)],'At Response', 'Pre Contrast Change');
% set(lg, 'fontsize', 24, 'box', 'on', 'edgecolor', 'w', 'location', 'southeast');

%ylims=[-.2 1.1];
%ylims=[4.75 7.09]; 
sp(12)=subplot(2,7,12); hold on; box off; 
for z=1:3
    pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),abs(nanmean(PreEvSlopeContra(:,z,subs),3)-nanmean(PreEvSlopeIpsi(:,z,subs),3)),'linestyle', '-', 'linewidth',3, 'color', colours_red{z+2});
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(PreEvSlopeContra(:,z,subs),3)*-1,'linestyle', '-', 'linewidth',3, 'color', colours_red{z+2});
    %pl2(z)=plot(nanmean(Plot_RTBins(:,z,subs),3),nanmean(PreEvSlopeIpsi(:,z,subs),3)*-1,'linestyle', '--', 'linewidth',3, 'color', colours_red{z+2});
    l(z)=line([delays(z) delays(z)],[ylims(1) ylims(2)], 'linestyle', '-', 'linewidth',2, 'color', colours_red{z+2});
    %ll(z)=line([0 0],[ylims], 'linestyle', '-', 'linewidth',3, 'color', colours_sessions{5});
    line([500 3600], [0 0],'linestyle', '-', 'linewidth',2, 'color', 'k');
end
 xlim([500 3100]); ylim(ylims);
set(gca,'fontsize', 18, 'xtick', [500 1000 1500 2000 2500 3000 3500], 'xticklabel', {'' '1000' '' '2000' '' '3000' '' }, ...;
    'ytick', [ -.002  -.001 0 .001 .002 .003], 'yticklabels', {'-.002' '-.001' '0' '.001' '.002' '.003'});
ylabel('Slope','fontsize',24); xlabel('RT (ms)', 'fontsize',24);

set(sp(6), 'position', [-.05 .5 .25 .25]);
set(sp(7), 'position', [.2 .4 .3 .58]);
set(sp(13), 'position', [.4 .85 .1 .13]);
set(sp(8), 'position', [.52 .4 .11 .58]);
set(sp(14), 'position', [.52 .85 .05 .13]);
set(sp(10), 'position', [.71 .69 .12 .24]);
set(sp(9), 'position', [.87 .69 .12 .24]);
set(sp(12), 'position', [.71 .4 .12 .24]);
set(sp(11), 'position', [.87 .4 .12 .24]);

% h=axes; set(h,'position', [0 .01 .99 .99], 'visible', 'off');% Create a new axis that is invisible. Otherwise all text will be positioned relevant to your most recent lot which you might not want.
% text('String', 'Motor Preparation','fontsize', 30, 'fontweight', 'bold', 'position',[.01, .94]);
% text('String', 'H','fontsize', 34, 'fontweight', 'bold', 'position',[.01, .78]);
% text('String', 'I','fontsize', 34, 'fontweight', 'bold', 'position',[.16, .97]);
% text('String', 'J', 'fontsize', 34, 'fontweight', 'bold','position',[.51, .97]);
% text('String', 'K','fontsize', 34, 'fontweight', 'bold','position',[.67, .94]);
% text('String','L','fontsize', 34, 'fontweight', 'bold','position',[.85, .94]);
% text('String', 'M','fontsize', 34, 'fontweight', 'bold','position',[.67, .65]);
% text('String','N','fontsize', 34, 'fontweight', 'bold','position',[.85, .65]);

end 





%% Statistics
subs=[2:5 7:18 20:22];
Inputdata=squeeze(nanmean((BetaContraBins(:,:,subs) - BetaIpsiBins(:,:,subs)),1));
[T.BetaLatEv.H1 T.BetaLatEv.P1  T.BetaLatEv.CI1 T.BetaLatEv.Stats1]=ttest(Inputdata(1,:),0, .05/3); 
[T.BetaLatEv.H2 T.BetaLatEv.P2 T.BetaLatEv.CI2  T.BetaLatEv.Stats2]=ttest(Inputdata(2,:),0, .05/3); 
[T.BetaLatEv.H3 T.BetaLatEv.P3 T.BetaLatEv.CI3  T.BetaLatEv.Stats3]=ttest(Inputdata(3,:),0, .05/3); 
InputdataContra=squeeze(nanmean(BetaContraBins(:,:,subs),1)); InputdataIpsi= squeeze(nanmean(BetaIpsiBins(:,:,subs),1));
[T.BetaEv.H1 T.BetaEv.P1  T.BetaEv.CI1 T.BetaEv.Stats1]=ttest(InputdataContra(1,:),InputdataIpsi(1,:), .05/3); 
[T.BetaEv.H2 T.BetaEv.P2 T.BetaEv.CI2  T.BetaEv.Stats2]=ttest(InputdataContra(2,:),InputdataIpsi(2,:), .05/3); 
[T.BetaEv.H3 T.BetaEv.P3 T.BetaEv.CI3  T.BetaEv.Stats3]=ttest(InputdataContra(3,:),InputdataIpsi(3,:), .05/3); 
figure; 
for z=1:3 
    subplot(1,3,z); hold on; hist(Inputdata(z,:));
end

 Inputdata=squeeze(nanmean((PreEvSlopeContra(:,:,subs) - PreEvSlopeIpsi(:,:,subs)),1));
 [T.BetaSlope.H1 T.BetaSlope.P1  T.BetaSlope.CI1 T.BetaSlope.Stats1]=ttest(Inputdata(1,:),0, .05/3); 
 [T.BetaSlope.H2 T.BetaSlope.P2 T.BetaSlope.CI2  T.BetaSlope.Stats2]=ttest(Inputdata(2,:),0, .05/3); 
 [T.BetaSlope.H3 T.BetaSlope.P3 T.BetaSlope.CI3  T.BetaSlope.Stats3]=ttest(Inputdata(3,:),0, .05/3); 
% InputdataContra=squeeze(nanmean(PreEvSlopeContra(:,:,subs),1));InputdataIpsi= squeeze(nanmean(PreEvSlopeIpsi(:,:,subs),1));
% [T.BetaSlope.H1 T.BetaSlope.P1  T.BetaSlope.CI1 T.BetaSlope.Stats1]=ttest(InputdataContra(1,:),InputdataIpsi(1,:), .05/3); 
% [T.BetaSlope.H2 T.BetaSlope.P2 T.BetaSlope.CI2  T.BetaSlope.Stats2]=ttest(InputdataContra(2,:),InputdataIpsi(2,:), .05/3); 
% [T.BetaSlope.H3 T.BetaSlope.P3 T.BetaSlope.CI3  T.BetaSlope.Stats3]=ttest(InputdataContra(3,:),InputdataIpsi(3,:), .05/3); 
for z=1:3 
    subplot(1,3,z); hold on; hist(Inputdata(z,:));
end


clear Inputdata
Inputdata(:,:,1,:)= BetaIpsiBins(:,:,subs); 
Inputdata(:,:,2,:)= BetaContraBins(:,:,subs); 
ANOVA.Beta_RT_delay = RMAOV(Inputdata(:,:,:,:));
%figure; hist(ANOVA.Beta_RT_delay.Normality.StdResiduals2);

clear Inputdata
Inputdata(:,:,:)= PreEvSlopeContra(:,:,subs)-PreEvSlopeIpsi(:,:,subs); 
ANOVA.BetaLatSlope_RT_delay = RMAOV(Inputdata(:,:,:,:));
%%figure; hist(ANOVA.BetaSlope_RT_delay.Normality.StdResiduals2);

clear Inputdata
Inputdata(:,:,1,:)= BetaIpsiBins(:,:,subs); 
Inputdata(:,:,2,:)= BetaContraBins(:,:,subs); 
ANOVA.BetaLat_RT_delay = RMAOV(Inputdata(:,:,2,:)-Inputdata(:,:,1,:));
%figure; hist(ANOVA.BetaLat_RT_delay.Normality.StdResiduals2);

clear Inputdata
Inputdata(:,:,1,:)= BetaRIpsiBins(:,:,subs); 
Inputdata(:,:,2,:)= BetaRContraBins(:,:,subs); 
ANOVA.Betar_RT_delay = RMAOV(Inputdata);
ANOVA.BetarContra_RT_delay = RMAOV(Inputdata(:,:,2,:));
clear Inputdata
Inputdata(:,:,1,:)= BetaRIpsiBins(:,:,subs); 
Inputdata(:,:,2,:)= BetaRContraBins(:,:,subs); 
ANOVA.BetarLat_RT_delay = RMAOV(Inputdata(:,:,2,:)-Inputdata(:,:,1,:));



clear Inputdata
Inputdata(:,:,1,:)= PreEvSlopeIpsi(:,:,subs); 
Inputdata(:,:,2,:)= PreEvSlopeContra(:,:,subs); 
ANOVA.BetaSlope_PreEv_delay = RMAOV(Inputdata(:,3,:,:));%-Inputdata(:,:,1,:));
%%figure; hist(ANOVA.BetaSlope_PreEv_delay.Normality.StdResiduals2);

clear Inputdata
Inputdata(:,:,1,:)= PreEvSlopeIpsi(:,:,subs); 
Inputdata(:,:,2,:)= PreEvSlopeContra(:,:,subs); 
ANOVA.BetaLatSlope_PreEv_delay = RMAOV(Inputdata(:,:,2,:)-Inputdata(:,:,1,:));%-Inputdata(:,:,1,:));
%%figure; hist(ANOVA.BetaSlope_PreEv_delay.Normality.StdResiduals2);


clear Inputdata
Inputdata(:,:,1,:)= PreRTSlopeIpsi(:,1:3,subs); 
Inputdata(:,:,2,:)= PreRTSlopeContra(:,1:3,subs); 
ANOVA.BetaLatSlope_PreRT_delay = RMAOV(Inputdata(:,:,1,:)-Inputdata(:,:,2,:));
%%figure; hist(ANOVA.BetaSlope_PreRT_delay.Normality.StdResiduals2);

