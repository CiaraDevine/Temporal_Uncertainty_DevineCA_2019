 % SAIIT with 2AFC - press left if left-tilted grating gets brighter than right, press right if right brighter

% Contents:

% 1) load monitor params, set task parameters

% 2) Get filename

% 3) set up triggers and open PTB window

% 4) make stimuli

% 5) make stimulus sequences, save as vectors of textures.
% 6) designate trigger codes

% 7) make randomized sequence of trial types

% 8) present instructions to subject

% 9) Start trials



clear all;
 
commandwindow;

% reset(RandStream.getDefaultStream,sum(100*clock));

Screen('Preference', 'SkipSyncTests', 2);

path = {'C:\Users\eeglab\Desktop\Practice_code\Mats\'};

rewards(1,150:2001)=80:-0.0216:40; %% 0 points for RT<150 then points drop from 80 to 40 across target duration

% ***************************************************** BASIC SET - UP 

try


dlg_title = 'SAIIT Practice Conditions';

par.runID = 'Default';

    prompt = {'Enter SUBJECT/RUN/TASK IDENTIFIER:'};
    def = {par.runID};
    answer = inputdlg(prompt,dlg_title,1,def);
    par.runID = answer{1};

    
    try load([par.runID '.mat']); 
    par.recordEEG=1; par.useEL = 0; par.blocknumber =blck+1;
    found =1;
    catch

    par.runID='Default'; par.recordEEG=1; par.useEL = 0; par.blocknumber =1;
    found =0;
    end

monitorwidth_cm = 32;   % monitor width in cm

dist_cm = 57;  % viewing distance in cm



whichScreen = 0;    % Which screen (in lab there's two, and the one in the booth is number 2 according to windows dsiplay settings)

[scresw, scresh]=Screen('WindowSize',whichScreen);  % Get screen resolution

scres = [scresw, scresh]



center = [scresw scresh]/2;     % useful to have the pixel coordinates of the very center of the screen (usually where you have someone fixate)



hz=Screen('FrameRate', whichScreen,1);



disp(['SCREEN RESOLUTION ' num2str(scresw) ' x ' num2str(scresh)])

disp(['MONITOR REFRESH RATE ' num2str(hz) ' Hz']);



cm2px = scresw/monitorwidth_cm;  % multiplication factor to convert cm to pixels

deg2px = dist_cm*cm2px*pi/180;      % multiplication factor to convert degrees to pixels (uses aproximation tanT ~= T).



load gammafn_TCD        % parameters for the monitor's gamma function (relates luminance to brightness level 0-255)
gam=2.65;
maxLum = GrayLevel2Lum(255,Cg,gam,b0);  % get maximum luminance. For making stimuli like gratings, we often define the stimulus in terms of luminance values and then use one of these translating functions to convert to 0-255 at the very end.



fixRect = [center-2 center+2];

par.BGcolor = Lum2GrayLevel(maxLum/2,Cg,gam,b0); % mid-gray background



%%%%%%%%% IMPORTANT SETTINGS (some of these are asked for in the dialog, so here they serve as defaults)
if found == 0;
par.MAXcontrast = 0.8;      % Asymptotic Contrast - it will eventually end up there
end

par.RampTime = 0.1; % in seconds



par.numtrials = 50;            % Number of targets

par.delays = [800 1200 1600];       % delays between cue onset and target onset in ms



% Other Settings

par.leadintime = 1000;              % how long to pause before experiment starts



% *******************************************************************************************

% ASK EXPERIMENTER TO INPUT SESSION/BLOCK (LOGFILE) NAME AND ADJUST EXPERIMENTAL PARAMETERS

% GUI will come up and experiment will start loading after you click 'OK' in GUI

dlg_title = 'SAIIT Practice Conditions';

while 1

    prompt = {'Enter SUBJECT/RUN/TASK IDENTIFIER:','EEG? (1=yes, 0=no)','Maximum contrast','ramp time','block number'};

    def = {par.runID,num2str(par.recordEEG),num2str(par.MAXcontrast),num2str(par.RampTime),num2str(par.blocknumber)};

    answer = inputdlg(prompt,dlg_title,1,def);

    par.runID = answer{1};

    par.recordEEG = str2num(answer{2});

    par.MAXcontrast = str2num(answer{3});

    par.RampTime = str2num(answer{4});
    
    par.blocknumber = str2num(answer{5});

  

    %if exist([par.runID '.mat'],'file'),

     %   dlg_title = [par.runID '.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH']

   % else

        break;

    %end

end



% **********************************************************************************
ioObj=io64;%%
 status=io64(ioObj);
% Define a whole lot of characteristics of the physical stimulus

par.videoFrate = 100;   % Monitor refresh rate

par.FlickF = [20 25];      % Flicker frequencies of two sequency assignment to left and right

% Note the left-tilted stimulus is always stimulus "1" and right-tilted stimulus "2"

par.spatfreq = 1;       % Spatial frequency of gratings

par.outerrad_deg = 6;   % in DEGREES

par.innerrad_deg = 1;   % in DEGREES

par.targdur = 2;  % in sec. Return ramp will be at double rate. Choose multiple of 0.4!

par.BLcontrast = 0.5;       % contrast at baseline



% Set up for triggers

if par.recordEEG

    % Parallel Port for triggers - set to zero to start

   port = hex2dec('EC00');  % TRINITY

%     port = hex2dec('1130');   % CCNY

    io64(ioObj,port,0);

end



if par.useEL, ELCalibrateDialog, end % Calibrate eyeLink



% Opens a graphics window on the main monitor

window = Screen('OpenWindow', whichScreen, par.BGcolor);



if scres ~=[1024,768]
   error(['The monitor is NOT SET to the Screen Resolution. Please change it.'])
end

% if abs(hz-par.videoFrate)>1
% 
%     error(['The monitor is NOT SET to the desired frame rate of ' num2str(par.videoFrate) ' Hz. Change it.'])
% 
% end



if par.useEL

    %%%%%%%%% EYETRACKING PARAMETERS

    par.FixWinSize = 3;    % RADIUS of fixation (circular) window in degrees

    par.TgWinSize = 3;    % RADIUS of fixation (circular) window in degrees

    ELsetupCalib

    Eyelink('Command', 'clear_screen 0')

%     Eyelink('command', 'draw_box %d %d %d %d 15', center(1)-deg2px*par.FixWinSize, center(2)-deg2px*par.FixWinSize, center(1)+deg2px*par.FixWinSize, center(2)+deg2px*par.FixWinSize);

end



%  **********************  MAKE STIMULI

par.oriL = 135 * pi/180;

par.oriR = 45 * pi/180;

Rout = round(deg2px*par.outerrad_deg);  % radii in pix

Rin = round(deg2px*par.innerrad_deg);

D=Rout*2+1;                             % full stimulus size "D"

% Make a sinusoidal grating filling the stimulus rectangle:

[x,y] = meshgrid([1:D]-(D+1)/2,[1:D]-(D+1)/2);

clear GL GR

GL{1} = sin(par.spatfreq/deg2px*2*pi*(x.*cos(par.oriL)+y.*sin(par.oriL)));    % range -1 to 1 (needs to be transformed to brightness scale)

GR{1}= sin(par.spatfreq/deg2px*2*pi*(x.*cos(par.oriR)+y.*sin(par.oriR)));    % range -1 to 1 (needs to be transformed to brightness scale)

GL{2} = sin(par.spatfreq/deg2px*2*pi*(x.*cos(par.oriL)+y.*sin(par.oriL))+pi/2);    % pi/2 phase shift, range -1 to 1 (needs to be transformed to brightness scale)

GR{2} = sin(par.spatfreq/deg2px*2*pi*(x.*cos(par.oriR)+y.*sin(par.oriR))+pi/2);    % pi/2 phase shift, range -1 to 1 (needs to be transformed to brightness scale)



midLum = GrayLevel2Lum(par.BGcolor,Cg,gam,b0);   % The very middle luminance on the monitor in cd/m^2

lumAmpl = floor(midLum);   % luminance amplitude (divergence from midLum) in cd/m^2

for g = 1:2 % for both spatial phases... (will alternate to stave off adaptation a little...)

    GL{g}(find(GL{g}>0))=lumAmpl; GL{g}(find(GL{g}<0))=-lumAmpl;    % convert sinusoidal luminance modulation of the spatial pattern to square wave

    GR{g}(find(GR{g}>0))=lumAmpl; GR{g}(find(GR{g}<0))=-lumAmpl;

    % Cut out the annulus shape:

    for j=1:D

        for k=1:D

            [th,r]=cart2pol(x(j,k),y(j,k)); % cartesian to polar

            if r < Rin | r > Rout

                GL{g}(j,k)= 0; GR{g}(j,k)= 0;

            end

        end

    end

end

% Now we'll make frame sequences, which comprise just a vector of multipliers for the pattern stimuli we've generated

% above (mostly 0 and 1 for off and on, and -1 for reversed pattern)

framesperflickercycle = round(par.videoFrate./par.FlickF);



BLframeseq = []; TGframeseq = [];



% baseline frame sequence length (this will be repeated again and again in the ITI)

BLframeseqlen = LCM_SK(framesperflickercycle)*2; % double it because every second "cycle" is phase reversed - for 20Hz and 25Hz stimulation on 100Hz monitor this makes a 40-frame or 400-ms complete single cycle of the stimulus

BLframeseqlen_ms = BLframeseqlen*1000/par.videoFrate;



numrefr_intro = 80; % Length of intro ramping in (multiple of 40!)

rampin = [1:numrefr_intro]/numrefr_intro; % Ramp before cue onset

INTROframeseq = [];



% A standard baseline frame sequence (just columns of 1 and -1):

clear BLframeseq INTROframeseq

for f=1:length(par.FlickF)

    BLframeseq(:,f) = repmat([ones(1,framesperflickercycle(f)) -1*ones(1,framesperflickercycle(f))],1,BLframeseqlen/(2*framesperflickercycle(f)))';

    INTROframeseq(:,f) = repmat([ones(1,framesperflickercycle(f)) -1*ones(1,framesperflickercycle(f))],1,numrefr_intro/(2*framesperflickercycle(f)))';

end

BLframeseq(:,:,2) = [BLframeseq(11:end,:);BLframeseq(1:10,:)];      % SK_NEW

INTROframeseq(:,:,2) = [INTROframeseq(11:end,:);INTROframeseq(1:10,:)]; % SK_NEW



% NOW MAKE TEXTURES - these are the actual stimuli that will be presented

% make the blended baseline frame sequence specifically for 2 gratings superimposed:

% There are two different stimulus types, which are just shifted versions of one another (to counteract adaptation a

% little)

bCIF=par.BLcontrast;

clear introstim

for g=1:2   % two stimulus spatial phases

    for tp=1:2  % two different temporal phases % SK_NEW

        for n=1:size(BLframeseq,1)

            clear stim

            stim = midLum + (bCIF)*BLframeseq(n,1,tp)*GL{g}+(1-bCIF)*BLframeseq(n,2,tp)*GR{g};   % SK_NEW just the 3rd dimension index

            % Fixation point

            stim(Rout:Rout+2,Rout:Rout+2)=GrayLevel2Lum(255,Cg,gam,b0);

            BLstim(g,tp,n) = Screen('MakeTexture', window, Lum2GrayLevel(stim,Cg,gam,b0));

        end



        % Create stimuli for fade in

        for n=1:size(rampin,2)

            clear stim

            stim = midLum + ((bCIF)*INTROframeseq(n,1,tp)*GL{g}+(1-bCIF)*INTROframeseq(n,2,tp)*GR{g}).*rampin(n);     % SK_NEW just the 3rd dimension index



            introstim(g,tp,n) = Screen('MakeTexture', window, Lum2GrayLevel(stim,Cg,gam,b0));

        end

    end

end



% make targets:

numrefr = round(par.targdur*par.videoFrate);  % Target frame sequence length, ?? just the ramp down

ramprefr = round(par.videoFrate*par.RampTime);

CIF = [par.BLcontrast + [1:ramprefr]*(par.MAXcontrast-par.BLcontrast)./ramprefr ones(1,numrefr-ramprefr)*par.MAXcontrast];  % SK_NEW



TGframeseqlen = length(CIF);

lengthOK = [BLframeseqlen:BLframeseqlen:400];

if size(find(TGframeseqlen == lengthOK),2)<1;

    error(['Wrong stimulus length!'])

end

% Start with just the sequence of 1s and -1s again:

clear TGframeseq

for f=1:length(par.FlickF)

    TGframeseq(:,f) = repmat([ones(1,framesperflickercycle(f)) -1*ones(1,framesperflickercycle(f))],1,TGframeseqlen/(2*framesperflickercycle(f)))';

end

TGframeseq(:,:,2) = [TGframeseq(11:end,:);TGframeseq(1:10,:)];  % SK_NEW



% BLENDING for targets and make textures (dealing in actual luminance until last step):

for g = 1:2 % each of the two spatial phases

    for tp=1:2  % two different temporal phases % SK_NEW

        for n=1:size(TGframeseq,1)

            clear stim

            % A Left-tilt target:

            stim = midLum + (CIF(n))*TGframeseq(n,1,tp)*GL{g}+(1-CIF(n))*TGframeseq(n,2,tp)*GR{g};  % SK_NEW just indices

            % Fixation point

            stim(Rout:Rout+2,Rout:Rout+2)=GrayLevel2Lum(255,Cg,gam,b0);

            targstim(g,tp,n,1) = Screen('MakeTexture', window, Lum2GrayLevel(stim,Cg,gam,b0));   % SK_NEW just indices

            % A Right-tilt target:

            stim = midLum + (1-CIF(n))*TGframeseq(n,1,tp)*GL{g}+(CIF(n))*TGframeseq(n,2,tp)*GR{g};   % SK_NEW just indices



            targstim(g,tp,n,2) = Screen('MakeTexture', window, Lum2GrayLevel(stim,Cg,gam,b0));   % SK_NEW just indices

        end

    end

end



% PTB stimulus drawing functions need the stimulus rectangle on the screen

% to know where to put it... Here is coordinates centered on stim location:

stimrect = round([-1 -1 1 1]*D/2);



par.stimTPShift = repmat([1 1 2 2],1,ceil(par.numtrials/4));

par.stimTPShift(par.numtrials+1:end) = [];



%  ************************************************* CODES AND TRIAL SEQUENCE

% trigger codes - can only use these 15: [1 4 5 8 9 12 13 16 17 20 21 24 25 28 29]

par.CD_RESP  = 1;

par.CD_FIXON = 4;

par.CD_TGOFF = 5;   % target off

par.CD_TG = [8 9];   % target   % one for each target type

par.CD_BUTTONS = [12 13];

par.CD_BEEP = 29;           

par.CD_REWARD = 16;       % Feedback: reward

par.CD_PUNISH = 17;       % Feedback: punishment

par.CD_STARTBEEP = 28;

par.CD_CONDCUE = 21;

par.CD_FADEIN = 24; 



temp = repmat([1 2],[1,ceil(par.numtrials/2)]);

trialLR = temp(:,randperm(size(temp,2)));  % shuffle



% Things that we'll save on a trial by trial basis

clear ITIstartT TargOnT RespLR RespT resp

numResp=1;



% **************************************************************** START TASK   ************************************

% Instructions:

Screen('DrawText', window, 'Fixate on the central dot.', 0.05*scres(1), 0.15*scres(2), 255);

Screen('DrawText', window, 'Press LEFT button with LEFT hand as soon as', 0.05*scres(1), 0.35*scres(2), 255);

Screen('DrawText', window, 'you see that the LEFT-tilted pattern is stronger.', 0.05*scres(1), 0.45*scres(2), 255);

Screen('DrawText', window, 'Press RIGHT button with RIGHT hand as soon as', 0.05*scres(1), 0.60*scres(2), 255);

Screen('DrawText', window, 'you see that the RIGHT-tilted pattern is stronger.', 0.05*scres(1), 0.70*scres(2), 255);

% Screen('DrawTexture', window, targstim(s,trialLR(n)), [], [center center] + stimrect);

Screen('DrawText', window, 'Press to begin', 0.05*scres(1), 0.85*scres(2), 255);



Screen('Flip', window);

[clicks,x,y,whichButton] = GetClicks(whichScreen,0); % Waits for the user to press a button before starting

if par.recordEEG, io64(ioObj,port,par.CD_RESP); WaitSecs(0.005); io64(ioObj,port,0); end

if par.useEL, Eyelink('Message', ['TASK_START']); end

RespT(1) = GetSecs;

RespLR(1) = whichButton;  if RespLR(numResp)==3, RespLR(numResp)=2; end  % The first response will be the one that sets the task going, after subject reads instructions

resp.time = zeros(par.numtrials,1); resp.LR = zeros(par.numtrials,1);



% ***************************************************************************************************************

%%%%%%%%%%%%%%%%%%%% START TRIALS *******************************************************************************

% Initial lead-in:

Screen('FillRect',window, 255, fixRect);

Screen('Flip', window);
HideCursor;
WaitSecs(par.leadintime/1000);



% Start Task:

portUP=0; lastTTL=0; ButtonDown=0; nPTBevent = 0; PTBeventT = [];



for n=1:par.numtrials
    
    
    clicked = 0;
    % For every trial the following are presented:

    % 1. 200ms of gray background

    % 2. 80 frames of a baseline stimulus that fades in

    % 3. Show baseline stimulus (0.5 contrast for both patterns)

    % 4. Present target

    % 5. Show reward for last trial
% double click to start trial
        clear buttons
    while 1
        [xblah,yblah,buttons] = GetMouse(whichScreen);
        if length(find(buttons))>=2, break; end
    end

    % First show standard during ITI - the baseline frame sequence

    % 1. Fixation point

    Screen('FillRect',window, 255, fixRect);

    nPTBevent = nPTBevent+1;    % increment the number of events

    PTBevent(nPTBevent) = par.CD_FIXON;    % record the event that just happened (using the arbitrarily defined codes above)

    [VBLTimestamp PTBeventT(nPTBevent)] = Screen('Flip', window); 

    if (par.recordEEG), io64(ioObj,port,par.CD_FIXON); portUP=1; lastTTL=GetSecs; end

    % wait 200 ms (Natalie had it, dunno why)

    tstart = GetSecs;

    while GetSecs-tstart < 0.2

        if par.recordEEG, if portUP & GetSecs-lastTTL>0.01, io64(ioObj,port,0); portUP=0; end, end

    end

    g = rem(n,2)+1;

    par.stimShift(n) = g;

    

    % 2. Fading in flicker and NO cue in form of coloured fixation point

    nPTBevent = nPTBevent+1;    % increment the number of events

    PTBevent(nPTBevent) = par.CD_FADEIN;    % record the event that just happened (using the arbitrarily defined codes above)

    [VBLTimestamp PTBeventT(nPTBevent)] = Screen('Flip', window);

    if (par.recordEEG), io64(ioObj,port,par.CD_FADEIN); portUP=1; lastTTL=GetSecs; end

    for s=1:size(introstim,2)

        if par.recordEEG, if portUP & GetSecs-lastTTL>0.01, io64(ioObj,port,0); portUP=0; end, end

        Screen('DrawTexture', window, introstim(g,par.stimTPShift(n),s), [], [center center] + stimrect);

        Screen('FillRect',window, 255, fixRect);

        Screen('Flip', window);

    end

    

    % 3. delay period

    onsetdelay = par.delays(ceil(rand*3)); % choose random delay

    for m=1 : round(onsetdelay/BLframeseqlen_ms)

        for s=1:BLframeseqlen

            if par.recordEEG, if portUP & GetSecs-lastTTL>0.01, io64(ioObj,port,0); portUP=0; end, end

            Screen('DrawTexture', window, BLstim(g,par.stimTPShift(n),s), [], [center center] + stimrect);

            Screen('FillRect',window, 255, fixRect);

            if m==1 & s==1

                % Triggers to signal upcoming of the cue

                nPTBevent = nPTBevent+1;    % increment the number of events

                PTBevent(nPTBevent) = par.CD_CONDCUE; % Used to be: PTBevent(nPTBevent) = par.CD_FIXON;    % record the event that just happened (using the arbitrarily defined codes above)

                [VBLTimestamp PTBeventT(nPTBevent)] = Screen('Flip', window);

               % if (par.recordEEG), io64(ioObj,port,par.CD_CONDCUE); portUP=1; lastTTL=GetSecs; end

               % if (par.useEL), Eyelink('Message', ['TRIAL' num2str(n) 'CUE' num2str(par.CD_CONDCUE)]); end

            else

                Screen('Flip', window);

            end
        
           [ex,wy,buttons] = GetMouse(whichScreen);    % Get response 1=stimulation computer, 2=EEG booth 
        if any(buttons) && clicked==0 && m*BLframeseqlen_ms>300;
            if buttons(1), clicked=2; if (par.recordEEG), io64(ioObj,port,par.CD_BUTTONS(1)); portUP=1; lastTTL=GetSecs; end
            elseif buttons(end) && clicked==0, clicked=2; if (par.recordEEG), io64(ioObj,port,par.CD_BUTTONS(2)); portUP=1; lastTTL=GetSecs;end
            end
        end
            
        end

    end

    

    % 4. Present target
    time = GetSecs;

    before_time = GetSecs; b=0;
    tic;
    for s=1:TGframeseqlen

        if par.recordEEG, if portUP & GetSecs-lastTTL>0.01, io64(ioObj,port,0); portUP=0; end, end

        Screen('DrawTexture', window, targstim(g,par.stimTPShift(n),s,trialLR(n)), [], [center center] + stimrect);

        Screen('FillRect',window, 255, fixRect);

        if s==1

            nPTBevent = nPTBevent+1;    % increment the number of events

            PTBevent(nPTBevent) = par.CD_TG(trialLR(n));    % record the event that just happened (using the arbitrarily defined codes above)

            [VBLTimestamp PTBeventT(nPTBevent)] = Screen('Flip', window);

            if par.recordEEG, io64(ioObj,port,par.CD_TG(trialLR(n))); portUP=1; lastTTL=GetSecs; end

            if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) 'TG' num2str(par.CD_TG(trialLR(n)))]); end

        else

            Screen('Flip', window);

        end

       clear buttons
   
        [ex,wy,buttons] = GetMouse(whichScreen);    % Get response 1=stimulation computer, 2=EEG booth 
        if any(buttons) && clicked==0;
            if buttons(1), resp.LR(n) = 1; if (par.recordEEG && b == 0), io64(ioObj,port,par.CD_BUTTONS(1)); portUP=1; lastTTL=GetSecs; end

                if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) 'BTN' num2str(par.CD_BUTTONS(1))]); end

                nPTBevent = nPTBevent+1;    % increment the number of events

                PTBevent(nPTBevent) = par.CD_BUTTONS(1);    % record the event that just happened (using the arbitrarily defined codes above)

                PTBeventT(nPTBevent) = GetSecs;
                
                clicked = 1;

            elseif buttons(end), resp.LR(n) = 2; if (par.recordEEG && b == 0), io64(ioObj,port,par.CD_BUTTONS(2)); portUP=1; lastTTL=GetSecs; end;

                if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) 'BTN' num2str(par.CD_BUTTONS(2))]); end

                nPTBevent = nPTBevent+1;    % increment the number of events

                PTBevent(nPTBevent) = par.CD_BUTTONS(2);    % record the event that just happened (using the arbitrarily defined codes above)

                PTBeventT(nPTBevent) = GetSecs;
                
                clicked = 1;

            end
 

  
            time = GetSecs;

            if b == 0

                b = 1;

                resp.time(n) = time-before_time;

            end

        end

    end
    toc;
    
    nPTBevent = nPTBevent+1;    % increment the number of events

    PTBevent(nPTBevent) = par.CD_TGOFF;    % record the event that just happened (using the arbitrarily defined codes above)

    [VBLTimestamp PTBeventT(nPTBevent)] = Screen('Flip', window);

    if par.recordEEG, io64(ioObj,port,par.CD_TGOFF); portUP=1; lastTTL=GetSecs; end

    if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) 'TGOFF' num2str(par.CD_TGOFF)]); end

    

    % 5. Gray background with reward - 600ms

    Screen('Flip', window);

    % wait 400 ms (Natalie had it, dunno why)

    tstart = GetSecs;

    while GetSecs-tstart < 0.4

        if par.recordEEG, if portUP & GetSecs-lastTTL>0.01, io64(ioObj,port,0); portUP=0; end, end

    end
    
    if clicked==2;
 Screen('DrawText', window, 'Clicked Too Soon', 0.40*scres(1), 0.5*scres(2), 255);
 Screen('Flip', window);
    elseif clicked==1 && resp.time(n)<=0.15
 Screen('DrawText', window, 'Clicked Too Soon', 0.45*scres(1), 0.5*scres(2), 255);
 Screen('Flip', window);
    elseif resp.LR(n)==trialLR(n) && resp.time(n)>0.15
 Screen('DrawText', window, 'Well Done', 0.45*scres(1), 0.5*scres(2), 255);
 Screen('Flip', window);
    elseif resp.LR(n)~=trialLR(n) && resp.time(n)>0.15
 Screen('DrawText', window, 'Error', 0.45*scres(1), 0.5*scres(2), 255);
 Screen('Flip', window);
        elseif resp.time(n)==0
 Screen('DrawText', window, 'Too late', 0.45*scres(1), 0.5*scres(2), 255);
 Screen('Flip', window);
    end
    
end



if par.useEL,

    Eyelink('StopRecording');

    Eyelink('CloseFile');

    ELdownloadDataFile

end

%cleanup
Screen('CloseAll')

% Save all variables needed for later analysis

RespT = resp.time; RespLR = resp.LR;

%%%% calculate accuracy and RT for this block

% sometimes there is an extra trial tagged on to trialLR, need to figure
% out why
if length(trialLR) > length(resp.LR)
    trialLR=trialLR(1:(end-1));
end
    
accuracy(par.blocknumber) =(sum(resp.LR==trialLR' & resp.time>0.15 & resp.time<2)/length(trialLR))*100;
RTs = resp.time(resp.LR==trialLR' & resp.time>0.15 & resp.time<2);
medianRT(par.blocknumber) = median(RTs(find(zscore(RTs)<3 & zscore(RTs)>-3)));

disp(['Accuracy this block =', num2str(accuracy(end))])
disp(['Reaction Time this block =', num2str(medianRT(end))])

if par.blocknumber>1
disp(['Accuracy LAST block =', num2str(accuracy(end-1))])
disp(['Reaction Time LAST block =', num2str(medianRT(end-1))])
RTdiff = (medianRT(end)- medianRT(end-1))*1000;
disp(['RT diff =',num2str(RTdiff)]);
end

score(par.blocknumber) = sum(rewards(round(resp.time(resp.LR==trialLR' & resp.time<=2 & resp.time>0.15)*1000)));

blck = par.blocknumber;
id = par.runID;
%tconstant(par.blocknumber) = par.timeConstant;

save([par.runID],'PTBeventT','nPTBevent','PTBevent','RespT','RespLR','trialLR','par','resp', 'blck', 'accuracy', 'medianRT', 'score')  % save([par.runID],'ITIstartT','TargOnT','RespT','RespLR','trialITI','trialLR','par','resp')

save parPRAC par


xaxis = [0:(par.blocknumber+1)];
figure
subplot(2,2,1)
bar(accuracy);
title('Accuracy');
subplot(2,2,2)
bar(medianRT);
title('Reaction Time');
subplot(2,2,3)
bar(score);
title('Score');


catch   % If something went awry, print last error and exit:

    ShowCursor

    ple

    sca
end
    
