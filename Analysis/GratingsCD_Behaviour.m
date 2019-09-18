%% Delays Paper Behavioural Plots & Stats
clear all; close all; clc;

path = 'C:\Users\devineca\Desktop\Dropbox\Ciara_GratingsDelays\PreProcessing\Events\';
subjID =  { 'CJG'  'CR'  'AR'  'KD'  'TB'  'AOB'  'PM'  'CMG'  'SR'  'AB' 'GK' 'ED' 'SC' 'CB' 'CE' 'SH' 'RMC' 'JH', 'JR', 'ROC' 'KM' 'ROS'}; % In order of when they were tested 
subjects=22;
Delays=[800 1200 1600]; 
filename = [path 'BehaviouralOutput']; 
load(filename);

colours_sessions = {[1, .65, .67] [.95, .35, .40] [1, 0, .06] [.5, 0, 0] [0, 0, 0]}; % Red theme colours
colours_blocks ={[.46, 1, .58] [.28, 1, .4] [.28, 1, .88] [.14, .97, 1] [0, .86, .94] [.08, .63, 1] [.02, .4, 1] [.28, 0, .8] [.01, 0, .4] [0, 0, 0]}; %green/blue theme colours
grey = [.5, .5, .5]; 
subs=2:22;
%% Calculate accuracy, RT etc for each subject, delay and target direction 
RT2=RT+(Delay(1,:)*(1000/512));


for s=1:length(subjID)
    for z=1:3
        for b=1:50
            if ismember(sum(SubID==s & Block(2,:)==b & ~isnan(Score)), [48:54]);
                try
                    BlockScore(b,z,s)=nansum(datasample(Score(SubID==s &  Delay(2,:)==z & Block(2,:)==b),15));
                    numtrials(b,s)=sum(SubID==s & Block(2,:)==b);
                catch
                    BlockScore(b,z,s)=nan;
                    numtrials(b,s)=sum(SubID==s & Block(2,:)==b);
                end
            else
                BlockScore(b,z,s)=nan;
                numtrials(b,s)=sum(SubID==s & Block(2,:)==b);
            end
        end
    end

    for z=1:3 % delays      
        for c=1:2 %left vs right target
            idx= SubID==s & Targ==c+7 & Delay(2,:)==z; 
        Acc(c,z,s)=(nansum(Accuracy==1 & idx)/nansum(idx))*100; 
        rt(c,z,s) = nanmean(RT(idx)); 
        rt2(c,z,s) = nanmean(RT2(idx)); 
        CVar(c,z,s) = nanstd(RT(idx))/rt(c,z,s);
        score_sum(c,z,s) = sum(datasample(Score(idx& ~isnan(Score)),250)); 
        score_av(c,z,s) = nanmean(Score(idx& ~isnan(Score))); 
        Misses(c,z,s) = (nansum(Accuracy==2  & idx)/nansum(SubID==s & Delay(2,:)==z))*100; 
        Prem(c,z,s) = (nansum(RT<150 & idx)/nansum(SubID==s & Delay(2,:)==z))*100;
        RatioMP(c,z,s) = Misses(c,z,s)/Prem(c,z,s);
        
        if z==3 
            FastRT=[150];TooFast=-600; 
            idxp=SubID==s & Delay(2,:)==3 & RT<FastRT & RT>TooFast;
            NumPrem(s)=sum(idxp);
            ChoiceBias(s)=sum(idxp & Resp==12)/sum(idxp & Resp==13);
            if NumPrem(s)<10
                ChoiceBias(s)=nan;
            end
            if ChoiceBias(s)<1
            ChoiceBiasPlot(s)=(sum(idxp & Resp==13)/sum(idxp & Resp==12))*-1;
            else
            ChoiceBiasPlot(s)=ChoiceBias(s);
            end
        end
        end
    end
end

Data(:,1)=Sessid;
Data(:,2)=Accuracy;
Data(:,3)=RT;
Data(:,4)=Targ;
Data(:,5)=SubID;
Data(:,6)=Delay(2,:)';
save([path 'Learning_Model_Data']);


%% Choice Bias Test
[HBias PBias CIBias TBias]=ttest(ChoiceBias(subs),1);

figure; hold on;
spb=barh(1:length(subs), ChoiceBiasPlot(subs)*-1, 'k'); 
set(gca, 'fontsize', 24, 'box','off', 'ytick', [2:2:length(subs)]); 
ylabel('Subject', 'fontsize', 30); xlabel('Choice Bias (Premature Trials)', 'fontsize', 30);
line([0 0], [0 23], 'linewidth', 3, 'color', 'k');
xlim([-19 10]); ylim([0 23]);
%text(.5, -2, 'Left');
% %% Two-Way ANOVAs - including Target direction & delay 
% 
% Inputdata= Acc; 
% WSFactors= {'Targ' 'Delay'};
% WSLevels= {{ 'Left', 'Right'}...
%     {'800ms'; '1200ms'; '1600ms'}};
% WithinModel = ['Targ*Delay'];
% BSDesign= 0; Contrasts=0; 
% ANOVA.Accuracy.DelayXTarg = RMAOV(Inputdata); 
% 
% %hist(ANOVA.Acc.DelayXTarg.Normality.StdResiduals)
% 
% 
% Inputdata= rt; 
% WSFactors= {'Targ' 'Delay'};
% WSLevels= {{ 'Left', 'Right'}...
%     {'800ms'; '1200ms'; '1600ms'}};
% WithinModel = ['Targ*Delay'];
% BSDesign= 0; Contrasts=0; 
% ANOVA.RT.DelayXTarg = RMAOV(Inputdata); 
% %hist(ANOVA.RT.DelayXTarg.Normality.StdResiduals)
% 
% Inputdata= rt2; 
% WSFactors= {'Targ' 'Delay'};
% WSLevels= {{ 'Left', 'Right'}...
%     {'800ms'; '1200ms'; '1600ms'}};
% WithinModel = ['Targ*Delay'];
% BSDesign= 0; Contrasts=0; 
% ANOVA.RT2.DelayXTarg = RMAOV(Inputdata); 
% %hist(ANOVA.RT.DelayXTarg.Normality.StdResiduals)
% 
% Inputdata= CVar; 
% WSFactors= {'Targ' 'Delay'};
% WSLevels= {{ 'Left', 'Right'}...
%     {'800ms'; '1200ms'; '1600ms'}};
% WithinModel = ['Targ*Delay'];
% BSDesign= 0; Contrasts=0; 
% ANOVA.RTCVar.DelayXTarg = RMAOV(Inputdata); 
% %hist(ANOVA.RTCVar.DelayXTarg.Normality.StdResiduals)
% 
% Inputdata= Misses; 
% WSFactors= {'Targ' 'Delay'};
% WSLevels= {{ 'Left', 'Right'}...
%     {'800ms'; '1200ms'; '1600ms'}};
% WithinModel = ['Targ*Delay'];
% BSDesign= 0; Contrasts=0; 
% ANOVA.Misses.DelayXTarg = RMAOV(Inputdata); 
% %hist(ANOVA.Misses.DelayXTarg.Normality.StdResiduals)
% 
% Inputdata= Prem; 
% WSFactors= {'Targ' 'Delay'};
% WSLevels= {{ 'Left', 'Right'}...
%     {'800ms'; '1200ms'; '1600ms'}};
% WithinModel = ['Targ*Delay'];
% BSDesign= 0; Contrasts=0; 
% ANOVA.Prem.DelayXTarg = RMAOV(Inputdata); 
% %hist(ANOVA.Misses.DelayXTarg.Normality.StdResiduals)

%% ONE-WAY ANOVAs 

Inputdata= squeeze(nanmean(Acc(:,:,subs))); 
WSFactors= {'Delay'};
WSLevels= {{'800ms'; '1200ms'; '1600ms'}};
WithinModel = ['Delay'];
BSDesign= 0; Contrasts=0; 
ANOVA.Accuracy.Delay = RMAOV(Inputdata); 
ANOVAContrasts.Accuracy.Delay = RMAOV(Inputdata, [],[],[],1); 
[T.Acc.H1 T.Acc.P1  T.Acc.CI1 T.Acc.Stats1]=ttest(nanmean(Acc(:,1,subs)), nanmean(Acc(:,2,subs)), .05/3); 
[T.Acc.H2 T.Acc.P2 T.Acc.CI2  T.Acc.Stats2]=ttest(nanmean(Acc(:,1,subs)), nanmean(Acc(:,3,subs)), .05/3); 
[T.Acc.H3 T.Acc.P3 T.Acc.CI3  T.Acc.Stats3]=ttest(nanmean(Acc(:,2,subs)), nanmean(Acc(:,3,subs)), .05/3); 


Inputdata= squeeze(nanmean(score_av(:,:,subs))); 
WSFactors= {'Delay'};
WSLevels= {{'800ms'; '1200ms'; '1600ms'}};
WithinModel = ['Delay'];
BSDesign= 0; Contrasts=0; 
ANOVA.Score.Delay = RMAOV(Inputdata); 
ANOVAContrasts.Score.Delay = RMAOV(Inputdata, [],[],[],1); 
[T.Score.H1 T.Score.P1  T.Score.CI1 T.Score.Stats1]=ttest(nanmean(score_av(:,1,subs)), nanmean(score_av(:,2,subs)), .05/3); 
[T.Score.H2 T.Score.P2 T.Score.CI2  T.Score.Stats2]=ttest(nanmean(score_av(:,1,subs)), nanmean(score_av(:,3,subs)), .05/3); 
[T.Score.H3 T.Score.P3 T.Score.CI3  T.Score.Stats3]=ttest(nanmean(score_av(:,2,subs)), nanmean(score_av(:,3,subs)), .05/3); 


Inputdata= squeeze(nanmean(rt(:,:,subs))); 
WSFactors= {'Delay'};
WSLevels= {{'800ms'; '1200ms'; '1600ms'}};
WithinModel = ['Delay'];
BSDesign= 0; Contrasts=0; 
ANOVA.RT.Delay = RMAOV(Inputdata);  
ANOVAContrasts.RT.Delay = RMAOV(Inputdata,[],[],[],1);
[T.RT.H1 T.RT.P1  T.RT.CI1 T.RT.Stats1]=ttest(nanmean(rt(:,1,subs)), nanmean(rt(:,2,subs)), .05/3); 
[T.RT.H2 T.RT.P2 T.RT.CI2  T.RT.Stats2]=ttest(nanmean(rt(:,1,subs)), nanmean(rt(:,3,subs)), .05/3); 
[T.RT.H3 T.RT.P3 T.RT.CI3  T.RT.Stats3]=ttest(nanmean(rt(:,2,subs)), nanmean(rt(:,3,subs)), .05/3); 


Inputdata= squeeze(rt(:,:,subs)); 
ANOVA.RTLR.Delay = RMAOV(Inputdata);  

Inputdata= squeeze(nanmean(Misses(:,:,subs))); 
WSFactors= {'Delay'};
WSLevels= {{'800ms'; '1200ms'; '1600ms'}};
WithinModel = ['Delay'];
BSDesign= 0; Contrasts=0; 
ANOVA.Misses.Delay = RMAOV(Inputdata); 
ANOVAContrasts.Misses.Delay = RMAOV(Inputdata, [],[],[],1); 
%hist(ANOVA.Misses.Delay.Normality.StdResiduals)
[T.Misses.H1 T.Misses.P1  T.Misses.CI1 T.Misses.Stats1]=ttest(nanmean(Misses(:,1,subs)), nanmean(Misses(:,2,subs)), .05/3); 
[T.Misses.H2 T.Misses.P2 T.Misses.CI2  T.Misses.Stats2]=ttest(nanmean(Misses(:,1,subs)), nanmean(Misses(:,3,subs)), .05/3); 
[T.Misses.H3 T.Misses.P3 T.Misses.CI3  T.Misses.Stats3]=ttest(nanmean(Misses(:,2,subs)), nanmean(Misses(:,3,subs)), .05/3); 


Inputdata= squeeze(nanmean(Prem(:,:,subs))); 
WSFactors= {'Delay'};
WSLevels= {{'800ms'; '1200ms'; '1600ms'}};
WithinModel = ['Delay'];
BSDesign= 0; Contrasts=0; 
ANOVA.Prem.Delay = RMAOV(Inputdata); 
ANOVAContrasts.Prem.Delay = RMAOV(Inputdata, [],[],[],1); 
[T.Prem.H1 T.Prem.P1  T.Prem.CI1 T.Prem.Stats1]=ttest(nanmean(Prem(:,1,subs)), nanmean(Prem(:,2,subs)), .05/3); 
[T.Prem.H2 T.Prem.P2 T.Prem.CI2  T.Prem.Stats2]=ttest(nanmean(Prem(:,1,subs)), nanmean(Prem(:,3,subs)), .05/3); 
[T.Prem.H3 T.Prem.P3 T.Prem.CI3  T.Prem.Stats3]=ttest(nanmean(Prem(:,2,subs)), nanmean(Prem(:,3,subs)), .05/3); 

%% TWO-WAY ANOVAS ACCURACY VS delay & RT Bin
for s=1:length(subjID)
    PEarly(1,s) = sum(RT(Resp==12 & RT<150 & SubID==s))/sum(RT(SubID==s & RT<150))*100;
    PEarly(2,s) = sum(RT(Resp==13 & RT<150 & SubID==s))/sum(RT(SubID==s & RT<150))*100;
    
    for z=1:3
        for c=1:2
            
            numbins = 8;
            RTLims(:,c,z,s) = round(prctile(RT2(SubID==s & Delay(2,:)==z & Targ==c+7), linspace(0, 100, numbins+1)));
            RTLims2(:,z,s) = round(prctile(RT2(SubID==s & Delay(2,:)==z), linspace(0, 100, numbins+1)));
            for b=1:numbins
                AccRTBins2(b,z,s) = nansum(Accuracy==1 & SubID==s & Delay(2,:)==z & RT2>RTLims2(b, z,s) & RT2<RTLims2(b+1,z,s))/sum(SubID==s & Delay(2,:)==z & RT2>RTLims2(b,z,s) & RT2<RTLims2(b+1,z,s));
                ScoreRTBins2(b,z,s)= nanmean(Score(SubID==s & Delay(2,:)==z & RT2>RTLims2(b, z,s) & RT2<RTLims2(b+1,z,s)));
                RTBins2(b,z,s) =round(nanmean(RT2(SubID==s & Delay(2,:)==z & RT2>RTLims2(b,z,s) & RT2<RTLims2(b+1,z,s))));   
                numtrialsRTBins(b,z,s)=nansum(SubID==s & Delay(2,:)==z & RT2>RTLims2(b, z,s) & RT2<RTLims2(b+1,z,s));

                AccRTBins(b,c,z,s) = sum(Accuracy==1 & SubID==s & Delay(2,:)==z & Targ==c+7 & RT2>RTLims(b,c, z,s) & RT2<RTLims(b+1,c,z,s))/sum(SubID==s & Delay(2,:)==z & Targ==c+7 & RT2>RTLims(b,c,z,s) & RT2<RTLims(b+1,c,z,s));
                RTBins(b,c,z,s) =round(nanmean(RT(SubID==s & Delay(2,:)==z & Targ==c+7 & RT2>RTLims(b,c,z,s) & RT2<RTLims(b+1,c,z,s))));   
            end
        end
    end
end



Inputdata= AccRTBins2(:,:,subs); 
ANOVA.Accuracy.RTDelay = RMAOV(Inputdata); 
ANOVAContrasts.RTBin.Short = RMAOV(Inputdata(:,1,:),[],[],[],1); 
ANOVAContrasts.RTBin.Int = RMAOV(Inputdata(:,2,:),[],[],[],1); 
ANOVAContrasts.RTBin.Long = RMAOV(Inputdata(:,3,:),[],[],[],1); 


% % Check that results do not change if analysis done separately for left and
% % right targets 
% Inputdata= squeeze(AccRTBins(:,1,:,:)); 
% ANOVA.Accuracy.RTDelayLTarg = RMAOV(Inputdata); 
% 
% Inputdata= squeeze(AccRTBins(:,2,:,:)); 
% ANOVA.Accuracy.RTDelayRTarg = RMAOV(Inputdata); 
% %hist(ANOVA.Accuracy.RTDelay.Normality.StdResiduals)


%% Figure 1  behaviour 

for ploton1=1
f=figure; hold on;
temp=[1200 1200 1200];


sp(1) = subplot(2,4,1);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(round(mean(mean(Acc(:,:,subs),1),3))); plSE=round(std(mean(Acc(:,:,subs),1),[],3))./sqrt(length(subs));
errorbar(pldat, plSE,'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=3; add2=3.5; add3=6; add4=6.5;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '**','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '***','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{'800' '1200' '1600'},'fontsize', 20);
ylim([min(pldat)-2*add4 max(pldat)+2*add4]); ylabel('Accuracy (%)', 'fontsize', 24); 


sp(3) = subplot(2,4,3);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(round(mean(mean(rt(:,:,subs),1),3))); plSE=round(std(mean(rt(:,:,subs),1),[],3))./sqrt(length(subs));
errorbar(pldat, plSE,  'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=60; add2=70; add3=160; add4=170;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '***','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{'800' '1200' '1600'},'fontsize', 20);
ylim([min(pldat)-add4 max(pldat)+2*add4]);  ylabel('RT (ms)', 'fontsize', 24); 


sp(4) = subplot(2,4,4);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(mean(mean(Misses(:,:,subs),1),3)); plSE=std(mean(Misses(:,:,subs),1),[],3)./sqrt(length(subs));
errorbar(pldat, plSE,  'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=1; add2=1.1; add3=2; add4=2.1;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '***','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{'800' '1200' '1600'},'ytick', [0:1:8], 'yticklabels', {'' '1' '' '3' '' '5' '' '7' ''},'fontsize', 20);
ylim([0 max(pldat)+2*add4]); ylabel('Misses (%)', 'fontsize', 24); 

sp(5) = subplot(2,4,5);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(mean(mean(Prem(:,:,subs),1),3)); plSE=std(mean(Prem(:,:,subs),1),[],3)./sqrt(length(subs));
errorbar(pldat, plSE,  'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=1; add2=1.1; add3=2; add4=2.1;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '**','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '***','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{'800' '1200' '1600'},'ytick', [0:1:8], 'yticklabels', {'' '1' '' '3' '' '5' '' '7' ''},'fontsize', 20);
ylim([0 max(pldat)+2*add4]); ylabel('Prem Resp (%)', 'fontsize', 24); 
xlabel('Foreperiod', 'fontsize', 26);

sp(6) = subplot(2,4,6);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(round(mean(mean(score_av(:,:,subs),1),3))); plSE=round(std(mean(score_av(:,:,subs),1),[],3))./sqrt(length(subs));
errorbar(squeeze(round(mean(mean(score_av,1),3))), round(std(mean(score_av,1),[],3))./sqrt(length(subs)),  'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=3; add2=3.5; add3=6; add4=6.5;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '**','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{'800' '1200' '1600'},'fontsize', 20);
ylim([min(pldat)-add4 59]); ylabel('Points Per Trial', 'fontsize', 24); 


sp(7)= subplot(2,4,7);  hold on; box off;
ymaxx=3000;
for z=1:3
    histbin=temp(z):100:(Delays(z)+2000);
    plot(histbin, hist(RT2(Delay(2,:)==z), histbin), 'linewidth',4, 'color', colours_sessions{z+2});
    line([Delays(z) Delays(z)], [0 ymaxx],  'linewidth',2, 'color', colours_sessions{z+2});
      set(get(gca,'Ylabel'), 'color', 'k','string', 'Trial Count', 'fontsize', 22);
    set(get(gca,'Xlabel'), 'color', 'k','string', 'RT (ms)', 'fontsize', 22);
set(gca, 'linewidth', 1.5,'xlim', [600 3400], 'ylim', [0 2800], 'ycolor', 'k', 'fontsize', 20, 'xtick', [600 800 1200 1600 2000 2400 2800 3200 3600],'xticklabels',{'0..' '..800' '1200' '1600' '2000' '2400' '2800' '3200' '3600'}, 'ytick', [500:500:ymaxx]);

end

sp(8)= subplot(2,4,8);  hold on; box off;
for z=1:3
    hold on;
    histbins=25;
    %set(gca, 'linewidth', 1.5,'position', [.1 .09 .8 .9]);
    %plot(linspace(min(RT(Delay(2,:)==z)+Delays(z)),2000+Delays(z),histbins), hist(round(RT(Delay(2,:)==z)+Delays(z)), histbins),'linestyle', '--', 'linewidth',2, 'color', colours_sessions{z+2});
     % [Ax H1(z) H2] = plotyy(nanmean(RTBins2(:,z,:),3)+Delays(z), nanmean(AccRTBins2(:,z,:),3)*100,linspace(min(RT(Delay(2,:)==z)+Delays(z)),2000+Delays(z),histbins), hist(round(RT(Delay(2,:)==z)+Delays(z)), histbins));
    H1(z)=errorbar(nanmean(RTBins2(:,z,:),3), nanmean(AccRTBins2(:,z,:),3)*100, nanstd((AccRTBins2(:,z,:)*100),[],3)/sqrt(length(subs)),'linestyle', '-', 'linewidth',4, 'color', colours_sessions{z+2});
    set(gca, 'linewidth', 1.5,'xlim', [600 3400], 'ylim', [39 99],'ycolor', 'k', 'fontsize', 20, 'xtick', [600 800 1200 1600 2000 2400 2800 3200 3600],'xticklabels',{}, 'ytick', [40:10:100]);
    %set(Ax(2),'ycolor', 'w', 'fontsize', 20,'xlim', [600 3600], 'ylim', [0 2550], 'xtick', [], 'xticklabels',{}, 'ytick', [500 1000 1500 2000 2500], 'yticklabels', []);
    set(get(gca,'Ylabel'), 'color', 'k','string', 'Accuracy (%)', 'fontsize', 22);
    %set(get(gca,'Xlabel'), 'color', 'k','string', 'RT (ms)', 'fontsize', 22);
    set(H1(z), 'linestyle', '-', 'linewidth',4, 'color', colours_sessions{z+2});  
    %set(H2, 'linestyle', '--', 'linewidth', 3, 'color', colours_sessions{z+2});
    l=line([Delays(z) Delays(z)], [36 99], 'color', colours_sessions{z+2}, 'linewidth', 2);
end
lg=legend([H1(1) H1(2) H1(3) ], '800ms FP ', '1200ms FP', '1600ms FP'); 
set(lg, 'fontsize', 23,'location', 'southeast', 'box', 'off');

set(sp(1), 'position', [.12 .6 .14 .34]);
set(sp(3), 'position', [.36 .6 .14 .34]);
set(sp(4), 'position', [.04 .13 .14 .34]);
set(sp(5), 'position', [.25 .13 .14 .34]);
set(sp(6), 'position', [.46 .13 .14 .34]);
set(sp(7), 'position', [.7 .1 .29 .35]);
set(sp(8), 'position', [.7 .55 .29 .35]);

h(1)=axes; set(h(1),'position', [.08 .96 .1 .1], 'visible', 'off');label=text(.01, .0, 'A'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(3)=axes; set(h(3),'position', [.31 .96 .1 .1], 'visible', 'off');label=text(.01, .0, 'B'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(4)=axes; set(h(4),'position', [.01 .49 .1 .1], 'visible', 'off');label=text(.01, .0, 'C'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(5)=axes; set(h(5),'position', [.22 .49 .1 .1], 'visible', 'off');label=text(.01, .0, 'D'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(6)=axes; set(h(6),'position', [.42 .49 .1 .1], 'visible', 'off');label=text(.01, .0, 'E'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(7)=axes; set(h(7),'position', [.65 .91 .1 .1], 'visible', 'off');label=text(.01, .0, 'F'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(8)=axes; set(h(8),'position', [.65 .46 .1 .1], 'visible', 'off');label=text(.01, .0, 'G'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot

end

%h(10)=axes; set(h(10),'position', [.0 .59 .6 .4], 'visible', 'off'); imshow(taskschematic, 'InitialMagnification', 'fit');


%% Figure 1 Alt behaviour 

for ploton1=1
f=figure; hold on;
temp=[1200 1200 1200];


sp(1) = subplot(2,4,1);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(round(mean(mean(Acc(:,:,subs),1),3))); plSE=round(std(mean(Acc(:,:,subs),1),[],3))./sqrt(length(subs));
errorbar(pldat, plSE,'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=3; add2=3.5; add3=6; add4=6.5;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '**','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '***','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{'800' '1200' '1600'},'fontsize', 20);
ylim([min(pldat)-2*add4 max(pldat)+2*add4]); ylabel('Accuracy (%)', 'fontsize', 24); 
xlabel('Foreperiod (ms)', 'fontsize', 24); 

sp(3) = subplot(2,4,3);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(round(mean(mean(rt(:,:,subs),1),3))); plSE=round(std(mean(rt(:,:,subs),1),[],3))./sqrt(length(subs));
errorbar(pldat, plSE,  'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=60; add2=70; add3=160; add4=170;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '***','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{},'fontsize', 20);
ylim([min(pldat)-add4 max(pldat)+2*add4]);  ylabel('RT (ms) ', 'fontsize', 24); 


sp(4) = subplot(2,4,4);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(mean(mean(Misses(:,:,subs),1),3)); plSE=std(mean(Misses(:,:,subs),1),[],3)./sqrt(length(subs));
errorbar(pldat, plSE,  'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=1; add2=1.1; add3=2; add4=2.1;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '***','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{},'ytick', [0:1:8], 'yticklabels', {'' '1' '' '3' '' '5' '' '7' ''},'fontsize', 20);
ylim([0 max(pldat)+2*add4]); ylabel('Misses (%)', 'fontsize', 24); 

sp(5) = subplot(2,4,5);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(mean(mean(Prem(:,:,subs),1),3)); plSE=std(mean(Prem(:,:,subs),1),[],3)./sqrt(length(subs));
errorbar(pldat, plSE,  'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=1; add2=1.1; add3=2; add4=2.1;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '**','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '***','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{},'ytick', [0:1:8], 'yticklabels', {'' '1' '' '3' '' '5' '' '7' ''},'fontsize', 20);
ylim([0 max(pldat)+2*add4]); ylabel('Premature Response (%)', 'fontsize', 24); 

sp(6) = subplot(2,4,6);  % sp is a 'handle' used below to reposition the plot
box 'off'; hold on;
pldat=squeeze(round(mean(mean(score_av(:,:,subs),1),3))); plSE=round(std(mean(score_av(:,:,subs),1),[],3))./sqrt(length(subs));
errorbar(squeeze(round(mean(mean(score_av,1),3))), round(std(mean(score_av,1),[],3))./sqrt(length(subs)),  'color','k', 'linestyle', '-', 'linewidth', 3, 'markerfacecolor', 'k', 'markersize', 2.5); % Linestyle none gets rid of the line connecting the error bars
add1=3; add2=3.5; add3=6; add4=6.5;
line([1.1 1.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([2.1 2.9], [max(pldat)+ add1 max(pldat)+ add1], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(2.35,max(pldat)+add2, '***','color', 'r', 'fontsize', 20);
line([1.1 2.9], [max(pldat)+ add3 max(pldat)+ add3], 'linestyle', '-', 'linewidth', 2, 'color', 'k'); text(1.9,max(pldat)+add4, '**','color', 'r', 'fontsize', 20);
xlim([.7 3.3]); set(gca, 'linewidth', 1.5,'xtick', [1:3], 'xticklabels',{'800' '1200' '1600'},'fontsize', 20);
ylim([min(pldat)-add4 59]); ylabel('Points Per Trial', 'fontsize', 24); 
xlabel('Foreperiod (ms)', 'fontsize', 24); 


sp(7)= subplot(2,4,7);  hold on; box off;
ymaxx=3000;
for z=1:3
    histbin=1000:100:(Delays(z)+2000);
    H(z)=plot(histbin, hist(RT2(Delay(2,:)==z & Accuracy==1), histbin), 'linewidth',4, 'color', colours_sessions{z+2});
    H1(z)=plot(histbin, hist(RT2(Delay(2,:)==z & Accuracy==1), histbin), 'linewidth',4, 'color', colours_sessions{z+2});
    Hh(z)=plot(histbin, hist(RT2(Delay(2,:)==z & ismember(Accuracy, [0 3])), histbin), 'linewidth',2, 'color', colours_sessions{z+2});
    line([Delays(z) Delays(z)], [0 ymaxx],  'linewidth',2, 'color', colours_sessions{z+2});
      set(get(gca,'Ylabel'), 'color', 'k','string', 'Trial Count', 'fontsize', 22);
    %set(get(gca,'Xlabel'), 'color', 'k','string', 'RT (ms)', 'fontsize', 22);
%set(gca, 'linewidth', 1.5,'xlim', [600 3400], 'ylim', [0 2800], 'ycolor', 'k', 'fontsize', 20, 'xtick', [600 800 1200 1600 2000 2400 2800 3200 3600],'xticklabels',{'0..' '..800' '1200' '1600' '2000' '2400' '2800' '3200' '3600'}, 'ytick', [500:500:ymaxx]);
set(gca, 'linewidth', 1.5,'xlim', [500 3599], 'ylim', [0 2150], 'ycolor', 'k', 'fontsize', 20, 'xtick', [500 800:400:3600],'xticklabels',{}, 'ytick', [500:500:ymaxx]);

end
lg=legend([H(1) H(2) H(3) H1(3) Hh(3)], '800ms FP ', '1200ms FP ', '1600ms FP ', 'Correct', 'Error'); 
set(lg, 'fontsize', 23,'location', 'southeast', 'box', 'off'); % , 'orientation', 'horizontal'


sp(8)= subplot(2,4,8);  hold on; box off;
for z=1:3
    hold on;
    histbins=25;
    %set(gca, 'linewidth', 1.5,'position', [.1 .09 .8 .9]);
    %plot(linspace(min(RT(Delay(2,:)==z)+Delays(z)),2000+Delays(z),histbins), hist(round(RT(Delay(2,:)==z)+Delays(z)), histbins),'linestyle', '--', 'linewidth',2, 'color', colours_sessions{z+2});
     % [Ax H1(z) H2] = plotyy(nanmean(RTBins2(:,z,:),3)+Delays(z), nanmean(AccRTBins2(:,z,:),3)*100,linspace(min(RT(Delay(2,:)==z)+Delays(z)),2000+Delays(z),histbins), hist(round(RT(Delay(2,:)==z)+Delays(z)), histbins));
    H1(z)=errorbar(nanmean(RTBins2(:,z,:),3), nanmean(AccRTBins2(:,z,:),3)*100, nanstd((AccRTBins2(:,z,:)*100),[],3)/sqrt(length(subs)),'linestyle', '-', 'linewidth',4, 'color', colours_sessions{z+2});
        %set(Ax(2),'ycolor', 'w', 'fontsize', 20,'xlim', [600 3600], 'ylim', [0 2550], 'xtick', [], 'xticklabels',{}, 'ytick', [500 1000 1500 2000 2500], 'yticklabels', []);
        set(gca, 'linewidth', 1.5,'xlim', [500 3599], 'ylim', [41 99],'ycolor', 'k', 'fontsize', 20, 'xtick', [500 800:400:3600],'xticklabels',...
        {'0..' '..800' '1200' '1600' '2000' '2400' '2800' '3200' '3600'},'ytick', [40:10:100]);
    set(get(gca,'Ylabel'), 'color', 'k','string', 'Accuracy (%)', 'fontsize', 22);
    set(get(gca,'Xlabel'), 'color', 'k','string', 'RT(ms)', 'fontsize', 24);
    set(H1(z), 'linestyle', '-', 'linewidth',4, 'color', colours_sessions{z+2});  
    %set(H2, 'linestyle', '--', 'linewidth', 3, 'color', colours_sessions{z+2});
    l=line([Delays(z) Delays(z)], [41 99], 'color', colours_sessions{z+2}, 'linewidth', 2);

end
lg=legend([H1(1) H1(2) H1(3) ], '800ms FP ', '1200ms FP ', '1600ms FP '); 
set(lg, 'fontsize', 23,'location', 'southeast', 'box', 'off', 'orientation', 'horizontal');
end


set(sp(1), 'position', [.52 .16 .18 .3]);
set(sp(3), 'position', [.52 .6 .18 .3]);
set(sp(4), 'position', [.818 .72 .18 .26]);
set(sp(5), 'position', [.818 .41 .18 .26]);
set(sp(6), 'position', [.818 .09 .18 .26]);
set(sp(7), 'position', [.075 .57 .32 .38]);
set(sp(8), 'position', [.075 .11 .32 .38]);

h(1)=axes; set(h(1),'position', [.0 .97 .1 .1], 'visible', 'off');label=text(.01, .0, 'B'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(3)=axes; set(h(3),'position', [.44 .92 .1 .1], 'visible', 'off');label=text(.01, .0, 'C'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(4)=axes; set(h(4),'position', [.0 .51 .1 .1], 'visible', 'off');label=text(.01, .0, 'D'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(5)=axes; set(h(5),'position', [.44 .49 .1 .1], 'visible', 'off');label=text(.01, .0, 'E'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(6)=axes; set(h(6),'position', [.735 .98 .1 .1], 'visible', 'off');label=text(.01, .0, 'F'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(7)=axes; set(h(7),'position', [.735 .68 .1 .1], 'visible', 'off');label=text(.01, .0, 'G'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot
h(8)=axes; set(h(8),'position', [.735 .36 .1 .1], 'visible', 'off');label=text(.01, .0, 'H'); set(label,'fontsize', 32, 'fontweight', 'bold'); % This line gives a label (A) to the current suplot


%% Premature responses 

[t.H t.P t.CI t.stats]=ttest(PEarly(1,:), PEarly(2,:));




% RT Distributons Left vs right 
leftrightcolours={[.2, .6, .5] [.91, .41, .17]};
leftrightcolours2={[.2, .1, .8] [.6, .2, .7]};

figure;  hold on; box off;
ymaxx=100;
for z=3
    histbins=1000:100:(Delays(z)+2000); 
    HLL=hist(RT2(ismember(Delay(2,:), z)  & Resp==12 & Targ==8), histbins);
    HLR=hist(RT2(ismember(Delay(2,:), z)  & Resp==12 & Targ==9), histbins);
    HRR=hist(RT2(ismember(Delay(2,:), z)  & Resp==13 & Targ==9), histbins);
    HRL=hist(RT2(ismember(Delay(2,:), z)  & Resp==13 & Targ==8), histbins);
    plot(histbins, (HLL./sum(Delay(2,:)==3))*100, 'linewidth',4, 'color', leftrightcolours{1});
    plot(histbins, (HLR./sum(Delay(2,:)==3))*100, 'linewidth',2, 'color', leftrightcolours{1});
    plot(histbins, (HRR./sum(Delay(2,:)==3))*100, 'linewidth',4, 'color', leftrightcolours{2});
    plot(histbins, (HRL./sum(Delay(2,:)==3))*100, 'linewidth',2, 'color', leftrightcolours{2});
    line([Delays(z) Delays(z)], [0 ymaxx],  'linewidth',2, 'color', colours_sessions{z+2});
    set(get(gca,'Ylabel'), 'color', 'k','string', 'Proportion (%)', 'fontsize', 30);
    set(get(gca,'Xlabel'), 'color', 'k','string', 'RT (ms)', 'fontsize', 30);
    set(gca, 'linewidth', 1.5,'xlim', [600 3400], 'ylim', [0 6], 'ycolor', 'k', 'fontsize', 24, 'xtick', [600 800 1200 1600 2000 2400 2800 3200 3600],'xticklabels',{'0..' '..800' '1200' '1600' '2000' '2400' '2800' '3200' '3600'}, 'ytick', [0:1:ymaxx]);
    lg=legend('Left Resp - Correct', 'Left Resp - Error', 'Right Resp - Correct', 'Right Resp - Error'); 
    set(lg, 'fontsize', 25, 'box', 'off');
end


figure; hold on; 
b=bar([1:3], [nanmean(rt(1,1,subs)) nanmean(rt(2,1,subs)); ...
    nanmean(rt(1,2,subs)) nanmean(rt(2,2,subs));...
    nanmean(rt(1,3,subs)) nanmean(rt(2,3,subs))], 'grouped');
b(1).FaceColor=leftrightcolours2{1};
b(2).FaceColor=leftrightcolours2{2};
ylabel('RT(ms)', 'fontsize', 30); 
xlabel('Foreperiod Duration', 'fontsize', 30); 
set(gca, 'fontsize', 24, 'xtick', [1 2 3], 'xticklabels', {'800ms' '1200ms' '1600ms'});
lg=legend('Left Target', 'Right Target'); set(lg, 'fontsize', 26, 'box', 'off');

figure; hold on; 
for t=1:2
    eb=errorbar([1:3], nanmean(rt(t,:,subs),3), nanstd(rt(t,:,subs),[],3)./sqrt(length(subs)));
    set(eb,  'linewidth', 3, 'color', leftrightcolours2{t});
end
ylabel('RT(ms)', 'fontsize', 30); 
xlabel('Foreperiod Duration', 'fontsize', 30); 
set(gca, 'fontsize', 24, 'xtick', [1 2 3], 'xticklabels', {'800ms' '1200ms' '1600ms'}, 'ytick',[600:200:1600]);
ylim([501 1599]);
lg=legend('Left Target', 'Right Target'); set(lg, 'fontsize', 26, 'box', 'off');


