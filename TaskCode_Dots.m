%% %%%%%%%%%%%%%%%%%%% Dots paradigm with stimulus - evidence onset delays  %%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%% BASIC SET - UP %%%%%%%%%%%%%%%%%
clear all;
commandwindow;
Screen('Preference', 'SkipSyncTests', 2);
rng('shuffle'); % reseed rand functions

try
    
    %% %%%%%%%%%%%%%%%% First Dialogue Box %%%%%%%%%%%%%%%
    dlg_title = 'Dot Motion Delays Task';
    par.runID = 'INITIALS_DelaysDots';
    prompt = {'Enter INITIALS/Task:'};
    def = {par.runID};
    answer = inputdlg(prompt,dlg_title,1,def);
    par.runID = answer{1};
    
    
    %% %%%%%%%%%%%%%% 2nd Dialogue Box %%%%%%%%%%%%%%%
    
    try load([par.runID '.mat']); block=block+1;
    catch; block=1; par.cohLevels =50;
    end
    
    
    dlg_title = 'Dot Motion Delays Task';
    while 1
        prompt = {'Enter SUBJECT/Task','EEG? (1=yes, 0=no)', 'Number of Trials (multiples of 6)', 'Block', 'Booth(1/2/3)', 'DotCoh'};
        def = {par.runID,num2str(1), num2str(54), num2str(block), '1', num2str(par.cohLevels)};
        answer = inputdlg(prompt,dlg_title,1,def);
        par.runID = answer{1};
        par.recordEEG = str2num(answer{2});
        par.numtrials=str2num(answer{3});
        par.block = str2num(answer{4});
        booth=str2num(answer{5});
        par.cohLevels =str2num(answer{6});
        break
    end
    
    %% Screen settings
    
    monitorwidth_cm_Booth1Booth3 = [32 32 32];   % Monitor width for each booth1 in eeglab
    monitorwidth_cm = monitorwidth_cm_Booth1Booth3(booth); % monitor width in cm
    dist_cm = 57;  % viewing distance in cm
    whichScreen = 0;    % Which screen (in lab there's two, and the one in the booth is number 2 according to windows dsiplay settings)
    [scresw, scresh]=Screen('WindowSize',whichScreen);  % Get screen resolution
    scres = [scresw, scresh];
    center = [scresw scresh]/2;     % useful to have the pixel coordinates of the very center of the screen (usually where you have someone fixate)
    hz=Screen('FrameRate', whichScreen,1);
    disp(['SCREEN RESOLUTION ' num2str(scresw) ' x ' num2str(scresh)])
    disp(['MONITOR REFRESH RATE ' num2str(hz) ' Hz']);
    cm2px = scresw/monitorwidth_cm;  % multiplication factor to convert cm to pixels
    deg2px = dist_cm*cm2px*pi/180;      % multiplication factor to convert degrees to pixels (uses aproximation tanT ~= T).
    fixRect = [center-2 center+2];
    
    par.videoFrate = 60;   % Monitor refresh rate
    par.FlickFdots = 60;      % Flicker frequency in Hz for the dots
    par.skipframes = 3;
    
    
    %% %%%%%%%% EXPERIMENT SETTINGS %%%%%%%%%%%%%
    
    par.numtrials = par.numtrials;
    par.cohLevels = par.cohLevels; % 0 10 20 40 coherence levels - these will be randomly interleaved - separate by spaces or commas [0 10 20 40] [0 80]
    par.BGcolor=0;      % BACKGROUND COLOR
    par.dotspeed = 6*par.skipframes;       % in degrees per second
    par.dotsize = 4;%4   % in pixels
    par.numdots = 100;%60  % this will be the number in a square; the ones outside the circle will be taken away
    par.dotpatchsize = 8;   % degrees of visual angle - diameter of dot patch
    par.patchloc = [0 0]; % patch location coordinates [x y] in degrees relative to center of screen
    par.motionDir = [180 0];    % in degrees relative to positive x-axis (0=rightward)
    par.rgbGRAY=[221 221 221];
    
    par.totalDotDur = [2800 3200 3600]; % in ms - how long are the dots up in total, between the incoherent prelude, cues, delay, etc
    par.MotionOnsetDelay = [800 1200 1600];  % in ms - when does the motion become coherent? Can be more than 1 and will be randomized. If using a cue, setting it to double the cue onset time will promote a rhythm
    par.cohMotionOffset = [par.totalDotDur];  % in ms - when does the motion become incoherent - set to totalDotDur if it will stay coherent til the end
    par.leadintime = 2000; % how long to pause before experiment starts
    
    CohLevs = [40 30 25 20 15 13 10 8 7 6 5 4 3 2 1];
    stair=[1 1 1];
    
    rewards(1,150:2001)=80:-0.0216:40; %% 0 points for RT<150 then points drop from 80 to 40 across target duration
    
    par.useEL=0;    % 1 to use EyeLink and track the eyes. 0 otherwise
    par.FixWinSize = 3;    % RADIUS of fixation (circular) window in degrees
    
    %%  Set up for triggers
    
    
    
    if par.recordEEG
        
        % Parallel Port for triggers - set to zero to start
        if booth==3,
            port = hex2dec('EC00');  % TRINITY - this port is specific to the 64 bit operating system in booth 3
            io64(ioObj,port,0);
            
        elseif booth==1
            port = hex2dec('2010');  % this port works for the 32 bit operating system in booth 1
            lptwrite(port,0);
        end
    end
    
    
    %%  ************************************************* CODES AND TRIAL SEQUENCE
    % trigger codes - can only use these 15 for brain products: [1 4 5 8 9 12 13 16 17 20 21 24 25 28 29]  (don't ask!)
    par.CD_RESP  = 1;
    par.CD_FIX_ON =4;
    par.CD_DOTS_ON = 24;
    par.CD_COHMOTION_ON = [8 9 12 13 16 17]; % leftward rightward motion
    par.CD_BUTTONS = [20 25];   % left and right mouse
    par.CD_COHMOTION_OFF = 5;
    
    par.FlickTdots = round(par.videoFrate/par.FlickFdots);      % In number of video frames
    nrefON = par.FlickTdots;                       % number of refreshes where the dots are ON
    nrefOFF = 0;                      % number of refreshes where the dots are OFF
    
    par.numFr = round((par.totalDotDur/1000)*par.FlickFdots);
    par.cohMotionOnFr = round(par.MotionOnsetDelay*par.FlickFdots/1000)+1;
    par.cohMotionOffFr = round(par.cohMotionOffset*par.FlickFdots/1000)+1;
    
    
    %%  Construct a cell array called "coh" which has a cell for the timecourse of coherence on each delay and motion direction
    
    cond = 0;
    clear conditiondescrip coh trigger triggerFr correctDir
    for d=1:length(par.MotionOnsetDelay)
        for m=1:length(par.motionDir)
            cond=cond+1;
            coh{cond}(1:length(par.motionDir), 1:par.numFr(d)) = 0;
            coh{cond}(m, par.cohMotionOnFr(d):par.cohMotionOffFr(d)-1) = par.cohLevels;
        end
    end
    
    
    temp = repmat(1:6,[1,9]);
    % Now RANDOMIZE the sequence of trials
    temp = temp(randperm(length(temp)));      % jumble temp
    
    
    %% ***************************************************** START TASK
    % Instructions:
    leftmargin = 0.02;
    window = Screen('OpenWindow', whichScreen,  par.BGcolor);
    Screen('Flip', window);
    
    Screen('DrawText', window, ['Pay attention to the dots moving randomly on screen.'], leftmargin*scres(1), 0.05*scres(2), 255);
    Screen('DrawText', window, ['After a period of time, a proportion of the the dots will'], leftmargin*scres(1), 0.1*scres(2), 255);
    Screen('DrawText', window, ['begin to move coherently in a leftward OR rightward direction.'], leftmargin*scres(1), 0.15*scres(2), 255);
    Screen('DrawText', window, ['The length of time between the onset of the random dot motion stimulus '], leftmargin*scres(1), 0.25*scres(2), 255);
    Screen('DrawText', window, ['and the onset of coherent motion varies unpredictably from trial to trial. '], leftmargin*scres(1), 0.3*scres(2), 255);
    
    Screen('DrawText', window, ['From the time coherent motion onsets you will have 2 seconds to respond.'], leftmargin*scres(1), 0.4*scres(2), 255);
    Screen('DrawText', window, ['If the dots are moving LEFTwards,'], leftmargin*scres(1), .45*scres(2), 255);
    Screen('DrawText', window, ['click the LEFT mouse button with your LEFT thumb. '], leftmargin*scres(1), .5*scres(2), 255);
    Screen('DrawText', window, ['If the dots are moving RIGHTwards,'], leftmargin*scres(1), .55*scres(2), 255);
    Screen('DrawText', window, ['click the RIGHT mouse button with your RIGHT thumb.'], leftmargin*scres(1), .6*scres(2), 255);
    Screen('DrawText', window, ['Correct responses are awarded 40 points'], leftmargin*scres(1), 0.7*scres(2), 255);
    Screen('DrawText', window, [' PLUS a 0-40 point speed bonus.'], leftmargin*scres(1), 0.75*scres(2), 255);
    Screen('DrawText', window, ['Incorrect responses are awarded 0 points.'], leftmargin*scres(1), 0.8*scres(2), 255);
    
    Screen('DrawText', window, 'Press any button to begin task.', leftmargin*scres(1), 0.9*scres(2), 255);
    Screen('Flip', window);
    HideCursor;
    
    % Things that we'll save on a trial by trial basis
    clear PTBtrigT PTBtrig ClickT Click RespLR perf button2modir
    RespT=[];
    nPTBtrig=0;
    numResp=1;
    
    % Waits for the user to press a button.
    [clicks,x,y,whichButton] = GetClicks(whichScreen,0);
    if par.recordEEG && booth==3, io64(ioObj,port, par.CD_RESP); %WaitSecs(0.005); io64(ioObj,port,0);
    elseif par.recordEEG && booth==1,  lptwrite(port, par.CD_RESP);end %WaitSecs(0.005); lptwrite(port,0);
    if par.useEL, Eyelink('Message', ['TASK_START']); end
    ClickT(1) = GetSecs;
    Click(1)=whichButton(1);    % The first response will be the one that sets the task going, after subject reads instructions
    
    %% %%%%%%%%%%%%%%%%%% START TRIALS
    portUP=0; lastTTL=0; ButtonDown=0; nPTBevent = 0; PTBeventT = [];
    
    % initial lead-in:
    Screen('FillRect',window, par.BGcolor); % screen blank
    Screen('Flip', window);
    WaitSecs(par.leadintime/2000);
    
    
    cc=1;
    for n=1:par.numtrials
   %% Record some details for later use   
        clicked=0;
        trialCond(block,n) = temp(block,n); % 1 = 800ms delay left, 2=800ms delay right, 3=1200ms delay left, 4=1200ms right, 5 = 1600ms delay left, 6=1600ms delay right
        if ismember(trialCond(block,n), [1 2]);
            DelayF=par.cohMotionOnFr(1);
            Delay(block,n) = 1;
        elseif ismember(trialCond(block,n), [3 4]);
            DelayF=par.cohMotionOnFr(2);
            Delay(block,n) = 2;        
        elseif ismember(trialCond(block,n), [5 6]);
            DelayF=par.cohMotionOnFr(3);
            Delay(block,n) = 3;            
        end
        
        resp.time(block,n)=0;
        resp.LR(block,n)=0;
        
        if ismember(trialCond(block,n), [1 3 5])
            trialLR(block,n)=1;
        elseif ismember(trialCond(block,n), [2 4 6])
            trialLR(block,n)=2;
        end
        
 
        %% MAKE THE DOTS FOR THIS TRIAL
        clear dots
        numFr = size(coh{trialCond(n)},2);

        % First generate dots at random locations on each frame
        for i=1:numFr
            for d=1:par.numdots
                dots(d,:,i) = [(rand-0.5)*par.dotpatchsize (rand-0.5)*par.dotpatchsize];
            end
        end
        
        
        % then add the coherence by selecting dots to move in certain direction relative to
        % previous frame. A different random set is selected each frame.
        for i=par.skipframes+1:numFr
            r = randperm(par.numdots);
            for m=1:length(par.motionDir)
                ncd = round(par.numdots*coh{trialCond(block,n)}(m,i)/100);
                randsel = r(1:ncd);
                % for the selected dots, move them in a particular direction
                dots(randsel,1,i) = dots(randsel,1,i-par.skipframes)+cos(par.motionDir(m)*pi/180)*par.dotspeed/par.FlickFdots;         % x-coordinate
                dots(randsel,2,i) = dots(randsel,2,i-par.skipframes)-sin(par.motionDir(m)*pi/180)*par.dotspeed/par.FlickFdots;         % y-coordinate
                r(1:ncd)=[];
            end
            % if it's gone off to the left, wrap it around to the far right
            dots(find(dots(:,1,i)<par.dotpatchsize/2),1,i) = dots(find(dots(:,1,i)<par.dotpatchsize/2),1,i)+par.dotpatchsize;
            % if it's gone off to the right, wrap it around to the far left
            dots(find(dots(:,1,i)>par.dotpatchsize/2),1,i) = dots(find(dots(:,1,i)>par.dotpatchsize/2),1,i)-par.dotpatchsize;
            % if it's gone off to the bottom, wrap it around to the top
            %         dots(find(dots(:,2,i)<par.dotpatchsize/2),1,i) = dots(find(dots(:,2,i)<par.dotpatchsize/2),1,i)+par.dotpatchsize;
            % if it's gone off to the top, wrap it around to the bottom
            %         dots(find(dots(:,2,i)>par.dotpatchsize/2),1,i) = dots(find(dots(:,2,i)>par.dotpatchsize/2),1,i)-par.dotpatchsize;
        end
        % Finally, go through the dots and get rid of the dots falling outside the
        % circle - put them off the screen.
        for i=1:numFr
            for d=1:par.numdots
                if sqrt(sum(dots(d,:,i).^2)) > par.dotpatchsize/2
                    dots(d,:,i) = 2*center/deg2px + 0.01;
                end
            end
        end
        savealldots{block}{n}=dots;
        %while GetSecs<ITIstart+par.ITI/1000; end
        
        %%  START TRIAL
        
        if n==(par.numtrials/2)+1
            Screen('FillRect',window,255, fixRect);
            Screen('Flip', window, [], 1);  % "Dont clear" = 1 means it won't clear the framebuffer after flip (will  keep fixation point)
            tbreakover = GetSecs;
            while GetSecs-tbreakover < 2 % Wait a second before restarting
            end
        else WaitSecs(1)
        end
        
        clear buttons;
        WaitSecs(1)
        
        Screen('FillRect',window,255, fixRect);
        Screen('Flip', window, [], 1);  % "Dont clear" = 1 means it won't clear the framebuffer after flip (will  keep fixation point)
        nPTBevent = nPTBevent+1;    % increment the number of events
        PTBevent(block,nPTBevent) = par.CD_FIX_ON;    % record the event that just happened (using the arbitrarily defined codes above)
        if par.recordEEG && booth==3, io64(ioObj,port, par.CD_FIX_ON); %WaitSecs(0.005); io64(ioObj,port,0);
        elseif par.recordEEG && booth==1,  lptwrite(port, par.CD_FIX_ON); %WaitSecs(0.005); lptwrite(port,0);
        end
        
        WaitSecs(.4);
        
        
        portUP=0; lastTTL=0; ButtonDown=0;
        
        % Present DOT MOTION
        for i=1:numFr
            clear buttons
            
            %% Dots On & Delay Period
            Screen('DrawDots', window, dots(:,:,i)'*deg2px, par.dotsize, par.rgbGRAY, center+par.patchloc.*[1 -1]);
            Screen('FillRect',window, 255, fixRect);
            Screen('Flip', window); % refresh
            
            if i==1
                nPTBevent = nPTBevent+1;    % increment the number of events
                PTBevent(nPTBevent) = par.CD_DOTS_ON;    % record the event that just happened (using the arbitrarily defined codes above)
                [VBLTimestamp PTBeventT(block,nPTBevent)] = Screen('Flip', window);
                if par.recordEEG && booth==3, io64(ioObj,port, par.CD_DOTS_ON); %WaitSecs(0.005); io64(ioObj,port,0);
                elseif par.recordEEG && booth==1, lptwrite(port, par.CD_DOTS_ON); %WaitSecs(0.005); lptwrite(port,0);
                end
            end
            if i<DelayF
                %Check for and record early responses
                [ex,wy,buttons] = GetMouse(whichScreen);    % Get response 1=stimulation computer, 2=EEG booth
                if any(buttons) && clicked==0 ;
                    if buttons(1), clicked=2;  resp.LR(block,n) = 3; resp.time(block,n)=NaN;
                        if par.recordEEG && booth==3, io64(ioObj,port, par.CD_BUTTONS(1)); %WaitSecs(0.005); io64(ioObj,port,0);
                        elseif par.recordEEG && booth==1, lptwrite(port, par.CD_BUTTONS(1)); %WaitSecs(0.005); lptwrite(port,0);
                        end
                        
                    elseif buttons(end) && clicked==0, clicked=2;
                        if par.recordEEG && booth==3, io64(ioObj,port, par.CD_BUTTONS(2)); %WaitSecs(0.005); io64(ioObj,port,0);
                        elseif par.recordEEG && booth==1, lptwrite(port, par.CD_BUTTONS(2)); %WaitSecs(0.005); lptwrite(port,0);
                        end
                        
                    end
                end
                
            elseif i>DelayF-1
                %% Coherent Motion
                if i == DelayF
                    targtime=GetSecs;
                    nPTBevent = nPTBevent+1;    % increment the number of events
                    PTBevent(nPTBevent) = par.CD_COHMOTION_ON(trialCond(block,n));    % record the event that just happened (using the arbitrarily defined codes above)
                    [VBLTimestamp PTBeventT(nPTBevent)] = Screen('Flip', window);
                    if (par.recordEEG) && booth==3, io64(ioObj,port,par.CD_COHMOTION_ON(trialCond(block,n)));
                    elseif (par.recordEEG) && booth==1, lptwrite(port,par.CD_COHMOTION_ON(trialCond(block,n)));
                    end
                end
                
                % Check for and record responses
                clear buttons
                [ex,wy,buttons] = GetMouse(whichScreen);    % Get response 1=stimulation computer, 2=EEG booth
                if any(buttons) && clicked==0;
                    if buttons(1), resp.LR(block,n) = 1;
                        if (par.recordEEG &&  booth==3), io64(ioObj,port,par.CD_BUTTONS(1));
                        elseif (par.recordEEG && booth==1), lptwrite(port,par.CD_BUTTONS(1));
                        end
                        if par.useEL, Eyelink('Message', ['TRIAL' num2str(block,n) 'BTN' num2str(par.CD_BUTTONS(1))]); end
                        
                        nPTBevent = nPTBevent+1;    % increment the number of events
                        PTBevent(nPTBevent) = par.CD_BUTTONS(1);    % record the event that just happened (using the arbitrarily defined codes above)
                        PTBeventT(nPTBevent) = GetSecs;
                        clicked = 1;
                        resp.time(block,n)= GetSecs-targtime;
                        
                    elseif buttons(end), resp.LR(block,n) = 2;
                        if (par.recordEEG && booth==3), io64(ioObj,port,par.CD_BUTTONS(2));
                        elseif (par.recordEEG &&  booth==1), lptwrite(port,par.CD_BUTTONS(2));
                        end
                        if par.useEL, Eyelink('Message', ['TRIAL' num2str(block,n) 'BTN' num2str(par.CD_BUTTONS(2))]); end
                        nPTBevent = nPTBevent+1;    % increment the number of events
                        PTBevent(nPTBevent) = par.CD_BUTTONS(2);    % record the event that just happened (using the arbitrarily defined codes above)
                        PTBeventT(nPTBevent) = GetSecs;
                        clicked = 1;
                        resp.time(block,n)= GetSecs-targtime;
                    end
                    
                end
                
            end
            
        end
        
        %% Trial Over
        nPTBevent = nPTBevent+1;    % increment the number of events
        PTBevent(nPTBevent) = par.CD_COHMOTION_OFF;    % record the event that just happened (using the arbitrarily defined codes above)
        [VBLTimestamp PTBeventT(nPTBevent)] = Screen('Flip', window);
        if (par.recordEEG) && booth==3, io64(ioObj,port,par.CD_COHMOTION_OFF);
        elseif (par.recordEEG) && booth==1, lptwrite(port,par.CD_COHMOTION_OFF);
        end
        
        WaitSecs(.5);
        
        %% Present Feedback
        tfeedback = GetSecs;
        if clicked==2;
            Screen('DrawText', window, 'Clicked Too Soon', 0.40*scres(1), 0.5*scres(2), 255);
            Screen('Flip', window);
        elseif clicked==1 && resp.time(block,n)<=0.15
            Screen('DrawText', window, 'Clicked Too Soon', 0.45*scres(1), 0.5*scres(2), 255);
            Screen('Flip', window);
        elseif resp.LR(block,n)==trialLR(block,n) && resp.time(block,n)>0.15
            Screen('DrawText', window, 'Well Done', 0.45*scres(1), 0.5*scres(2), 255);
            Screen('Flip', window);
        elseif resp.LR(block,n)~=trialLR(block,n) && resp.time(block,n)>0.15
            Screen('DrawText', window, 'Error', 0.45*scres(1), 0.5*scres(2), 255);
            Screen('Flip', window);
        elseif resp.time(block,n)==0
            Screen('DrawText', window, 'Too late', 0.45*scres(1), 0.5*scres(2), 255);
            Screen('Flip', window);
        end
        
        while GetSecs-tfeedback < 1 %% Leave feedback on for 1 second
        end
        
        
        %% Give a break mid way through each block
        
        if n==par.numtrials/2;
            Screen('DrawText', window, 'Take a Break.', 0.45*scres(1), 0.4*scres(2), 255);
            Screen('DrawText', window, 'When Ready Click Left & Right Buttons Simultaneously To Continue.', 0.3*scres(1), 0.5*scres(2), 255);
            Screen('Flip', window);
            
            while 1
                [xblah,yblah,buttons] = GetMouse(whichScreen);
                if length(find(buttons))>=2,
                    break;
                end
            end
            
            Screen('FillRect',window,0, fixRect);
            Screen('Flip', window, [], 1);  % "Dont clear" = 1 means it won't clear the framebuffer after flip (will  keep fixation point)
            
            tbreakover = GetSecs;
            while GetSecs-tbreakover < 1.5 % Wait a second before restarting
            end
            
            
        end
        
        
        %% Give block feedback
        
        if n==par.numtrials;
            
            clear RTs
            par.Accuracy(block) =sum(resp.LR(block,:)==trialLR(block,:) & resp.time(block,:)>0.15 & resp.time(block,:)<2.1)/length(trialLR(block,:))*100;
            RTs = resp.time(resp.LR==trialLR(block,:) & resp.time(block,:)>0.15 & resp.time(block,:)<2.1);
            par.medianRT(block) = median(round(RTs(find(zscore(RTs)<3 & zscore(RTs)>-3))*1000));
            par.score(block) = sum(rewards(round(resp.time(resp.LR(block,:)==trialLR(block,:) & resp.time(block,:)<=2 & resp.time(block,:)>0.15)*1000)));
            
            
            Screen('DrawText', window, 'This Block is Over: Please Call the Experimenter.', 0.35*scres(1), 0.25*scres(2), 255);
            Screen('DrawText', window, 'Your Performance on This Block: ', 0.4*scres(1), 0.45*scres(2), 255);
            Screen('DrawText', window, ['Your Accuracy Was:   ' num2str(par.Accuracy(block)) '%'], 0.4*scres(1), 0.55*scres(2), 255);
            Screen('DrawText', window, ['And Your Average Reaction Time Was:   ' num2str(par.medianRT(block)) 'milliseconds'], 0.4*scres(1), 0.65*scres(2), 255);
            Screen('DrawText', window, ['Giving You a Score of:   ' num2str(round(par.score(block))) '  out of  ' num2str(sum(par.numtrials)*80) ' possible points'], 0.4*scres(1), 0.75*scres(2), 255);
            Screen('Flip', window);
            
            
            
            clear buttons
            while 1
                [xblah,yblah,buttons] = GetMouse(whichScreen);
                if length(find(buttons))>=.4,
                    break;
                end
            end
            
            
            %% Save all variables needed for later analysis
            
            save([par.runID],'PTBeventT','nPTBevent','PTBevent','resp', 'trialLR','par', 'savealldots', 'block', 'Delay')  % save([par.runID],'ITIstartT','TargOnT','RespT','RespLR','trialITI','trialLR','par','resp')
            
            if par.useEL,
                Eyelink('StopRecording');
                Eyelink('CloseFile');
                ELdownloadDataFile
            end
            
        end
        
    end
    
    
    %% cleanup
    Screen('CloseAll')
    
    
catch
    ShowCursor
    ple
    sca
end

