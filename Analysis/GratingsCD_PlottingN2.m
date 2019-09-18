%% %%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
load('C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\Analyses\Temporal_Uncertainty_Delays\PlottingData\N2_PlottingData'); 
delays=[800 1200 1600];
chanlocs = readlocs('cap128.loc');
fs=512;
numsubs=22;
subs=[2:5 7:18 20:22];


%% Make 6hz filter for visualising CPP
[B,A]=butter(4,(6*2)/fs); % creating filter for ERP waveforms


%% Waveforms N2
% RT Binned  N2
clear Plot* RTBins*
for s=1:numsubs
    if s~=19
        for z=1:3
            PlotN2(:,z,s)=squeeze(nanmean(N2(:, SubID==s & ismember(Delay(:,1)', z)),2));
            N2Amp(z,s)=nanmean(nanmean(N2(128:150, SubID==s & ismember(Delay(:,1)', z)),2),1);
        end
    else
    end
end

%% 
 
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


