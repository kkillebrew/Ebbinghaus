%1024X768 85Hz
%rawdata:
%1: Condition(1=Pursuit/Retinal Eccentricity/Box Sise Relative to Bar
%             2=E/RS 3=P/RS, 4=P/E, 5-8 held for other conditions; 9=Quicktime)
%2: Staircase (1=start small 2=start big)
%3: Iteration (Which particular one?)
%4: GrowthRate
%5: Answer (1=shrink, 2=grow);
%6: Accuracy (1=right 0=wrong)
%7: Reversal? (1=yes 0=no);
%8: Successful? (1=yes 0=no);

clear
Screen('Preference', 'SkipSyncTests', 1);
datafile = input('Enter Subject Code:','s');
datafile_full = sprintf('%s_full',datafile);
ListenChar(2);                                  %surpress keyboard output to command window and script
HideCursor;

eyetracking = 0;
checktiming = 0; % boolean, should we print info about the duration between Flips during animation?

a=KbName('a');                                  %define response keys
k=KbName('k');

if eyetracking==1
    ppd=32.62; %with 50cm 1024x768 at 73cm
else
    ppd=40;     %imac                                    %set ppd according to monitor size and viewer distance
end


% colors [R G B] 0-255
bgcolor = [0 0 0];       % background
scolor  = [255 255 255]; % square
rcolor  = [0 0 0];       % rectangle
fcolor  = [255 0 0];     % fixation spot
tcolor  = [255 255 255]; % text

sqtype = 'FillRect'; % the Screen command for drawing the context square: 'FillRect', 'FillOval', 'FrameRect', etc.
bartype = 'FillRect'; % the Screen command for drawing the target bar: 'FillRect', 'FillOval', 'FrameRect', 'DrawTexture', etc.

if eyetracking==1
    [w,rect]=Screen('OpenWindow',1,bgcolor);  %Use with 1024x768 Monitor
else
    [w,rect]=Screen('OpenWindow',0,bgcolor,[0,0,round(1024/32.62)*ppd,round(768/32.62)*ppd]);  %Use with anything else
end
xc=rect(3)/2;                                   %define screen center
yc=rect(4)/2;
Priority(9);                                    %allocate max memory

[nums, names]=GetKeyboardIndices;               %find keyboards (may have to change the numbers depending on setup)
con_ID=nums(1);

if eyetracking==1
    dev_ID =nums(2);
else
    dev_ID =nums(1);
end

if strcmp(bartype,'DrawTexture')
    % Enable alpha blending for contrast manipulations
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    imgpix = imread('roundedbar.tif'); % 128x128x2 (gray + alpha)
    imgalpha = imgpix(:,:,2);
    imgpix = imgpix(:,:,1)./255;
    imgpix = cat(3,imgpix*rcolor(1),imgpix*rcolor(2),imgpix*rcolor(3),imgalpha); % RGBA
    bartexid = Screen('MakeTexture',w,imgpix,[],1);
end

refresh=85;                                     %set to monitor refresh rate
duration=.4;                                    %animation duration
duration=duration/(1/refresh);
fixduration=.5;                                    %pre-animation fixation
fixduration=fixduration/(1/refresh);
fixpause=.5;
sgrowth= 3;                                      %rate of square growth/stimulus movement
smin=ppd*2;                                     %sets starting square size
rmin= 0.80 * smin;                                    %rectangle starting size
fr=ppd/10;                                      %fixation dot radius
xoffset=30;%xc/2;                                   %distance from center everything starts
yoffset=sgrowth*duration;                       %how far up some things start

fixjitterfactor=1*ppd;                          %Set randomization parameters for start locations
squarejitterfactor=.02;
rectjitterfactor=.02;

stair_repetitions=3;                                  %repetitions of each unique staircase.
trials_per_stair=25;%50;                                  %trials in each staircase.

stairtype_list =[1 2];
nstairs=length(stairtype_list);
iter_list=1:stair_repetitions;
niter=length(iter_list);
condition_list = [0 22.5 45 67.5 90];%[-9 10 11 12 13];%[10 11 12 13 9];%[1 2 3 4 9];                        % list of conditions to include
nconds=length(condition_list);

%break_trials = .25:.25:.75; % list of proportion of total trials at which to offer subject a self-timed break
break_trials = .1:.1:.9;    % list of proportion of total trials at which to offer subject a self-timed break

%values to step through in stair case
growth_list = linspace(-1.2,1.2,18);
%[-2 -1.5 -1.25 -1 -.75 -.55 -.35 -.2 -.05 .05 .2 .35 .55 .75 1 1.25 1.5 2];   

sc_list=[];
prevan_list=1;

for i=1:nconds
    for j=1:nstairs
        for l=1:niter
            if j==2 % changed from l==2 so that stair_repetitions=1 works
                sc_list(i,j,l)=length(growth_list);
            else
                sc_list(i,j,l)=1;
            end
            prevan_list(i,j,l)=1;
            place_list(i,j,l)=1;
        end
    end
end
% Sets up lists to keep track of current position in growth list for each
% staircase, as well as whether the previous answer in this staircase was
% correct or incorrect, as well as where to put a value in the reversal
% list.

rawdata=[];                                 %set up data files to be produced
Pos_v_Trace=[];

reversal_list=[];                           %lists all reversal points for each staircase

variables_list=repmat(fullyfact([nconds nstairs niter]),[trials_per_stair 1]);      %variable combinations for each trial
variables_list=sortrows(variables_list);

trial_list = 1:size(variables_list,1);                    %sets lists for trial reuse when fixation is broken.

com_list=[];                                              %count up succesfully completed trials


r1=3;                                                   %dot sizes for eye tracker calibration
r2=10;
acceptable_error=1.5*ppd;                               %acceptable distance from fixation eyes may be
acceptable_nofix=10;                                     %number of times fixation may be seen as broken before a trial is botched


if eyetracking==1;
    ivx=iViewXInitDefaults;
    Screen('FillRect',w,[0 0 128],rect);
    [results]=EyeTrackerCalibrate(ivx,acceptable_error,dev_ID,con_ID,w);   %run calib which runs valid, which runs calib if it fails
    
    Screen('FillRect',w,[0 0 0],rect);
end

Screen('TextSize',w,24);                                               %Directions
text='Fixate on the red dot at all times during the experiment.';
width=RectWidth(Screen('TextBounds',w,text));
Screen('DrawText',w,text,xc-width/2,yc-100,tcolor);
text='At the end of each trial,';
width=RectWidth(Screen('TextBounds',w,text));
Screen('DrawText',w,text,xc-width/2,yc-50,tcolor);
text='indicate if the black bar appeared to grow or shrink.';
width=RectWidth(Screen('TextBounds',w,text));
Screen('DrawText',w,text,xc-width/2,yc,tcolor);
text='Press "A" for shrink or "K" for Grow';
width=RectWidth(Screen('TextBounds',w,text));
Screen('DrawText',w,text,xc-width/2,yc+50,tcolor);
text='when the stimulus disappears.';
width=RectWidth(Screen('TextBounds',w,text));
Screen('DrawText',w,text,xc-width/2,yc+100,tcolor);
text='Press any key to begin.';
width=RectWidth(Screen('TextBounds',w,text));
Screen('DrawText',w,text,xc-width/2,yc+150,tcolor);

Screen('Flip',w);
WaitSecs(1);
KbWait(dev_ID);

if eyetracking==1
    [result, ivx]=iViewXComm('open', ivx);              %open connection
end

z=0;

while length(com_list)<length(variables_list)           %keep going till all trials are successfully completed
    z=z+1;                                              %count up trial number
    botch=0;                                                    %Resets trial success/failure
    
    thistrial = Sample(trial_list(trial_list~=0));                  %grab an uncompleted trial
    
    cidx = variables_list(thistrial,1);
    condition = condition_list(cidx);                     %Choose variable levels
    rawdata(z,1)=condition;
    sidx = variables_list(thistrial,2);
    stair=stairtype_list(sidx);
    rawdata(z,2)=stair;
    iidx= variables_list(thistrial,3);
    iteration=iter_list(iidx);
    rawdata(z,3)=iteration;
    rgrowth=growth_list(sc_list(cidx,sidx,iidx));       %set growth level based on position in that particular staircase
    rawdata(z,4)=rgrowth;
    
    fixjitterx=((-1)^randi(2))*rand*fixjitterfactor;            %set random start locations
    fixjittery=((-1)^randi(2))*rand*fixjitterfactor;
    squarejitter=(squarejitterfactor*rand)*smin;
    rectjitter=(rectjitterfactor*rand)*smin;
    
    % give subject break at certain trials...
    this_b = 0;
    for b = break_trials
        if z == round(b*length(variables_list))
            this_b = b;
            break
        end
    end
    if this_b
        % display break message
        text=sprintf('You have completed %d%% of the trials.',round(b*100));
        width=RectWidth(Screen('TextBounds',w,text));
        Screen('DrawText',w,text,xc-width/2,yc,tcolor);
        text='Press any key when you are ready to continue.';
        width=RectWidth(Screen('TextBounds',w,text));
        Screen('DrawText',w,text,xc-width/2,yc+50,tcolor);
        Screen('Flip',w);
        KbWait(dev_ID);
    end
    
    sx=xc-xoffset+fixjitterx;                                               %Set fixation starting location
    sy=yc+fixjittery;
    
    sx1=sx;                                                                 %set initial fixation values
    sy1=sy;
    
    Screen('FillOval',w,[255 0 0],...                           %draw fixation
        [sx-fr,sy-fr,sx+fr,sy+fr]);
    
    xpoint=sx;                                                  %set starting points for eyetracker gaze checking
    ypoint=sy;
    xeye_cur = xpoint;
    yeye_cur = ypoint;
    
    Screen('Flip',w);
    
    if eyetracking==1
        fixcount=0;
        result=iViewXComm('send', ivx, 'ET_STR');           %datastreamingon
        while fixcount<fixduration/(1/85)                                       %set required continuous fixation time to start
            
            [xeye_new, yeye_new]=EyeTrackerGazeCheck(ivx,xeye_cur,yeye_cur,w);  %Track the eye position
            xeye_cur = xeye_new;
            yeye_cur = yeye_new;
            xdiff = xpoint-xeye_new;
            ydiff = ypoint-yeye_new;
            delta=sqrt((xdiff^2)+(ydiff^2));                                    %get the difference in pixels between the fixation and eye position
            
            if delta>acceptable_error
                fixcount=0;                                                    %reset every time fixation is broken
            else
                fixcount=fixcount+1;                                                    %count number of continuous successful fixations
            end
            
        end
        result=iViewXComm('send', ivx, 'ET_EST'); %datastreamingoff
    else
        WaitSecs(.5);  %when not using eyetracking
    end
    
    switch condition                                                    %set initial square and bar values
        case {1 2 3 9}
            % central fixation, small square centered, bar left of fixation
            srect = [sx1-smin sy1-smin sx1+smin sy1+smin];
            srect = srect + [-squarejitter -squarejitter squarejitter squarejitter]; % add size jitter
            rrect = [sx1-smin+(ppd*.25) sy1-rmin sx1-smin+(ppd*.75) sy1+rmin];
            rrect = rrect + [-rectjitter -rectjitter rectjitter rectjitter]; % add size jitter
        case {10 11 12 13 -9 0 22.5 45 67.5 90}
            % central fixation, small square centered, bar right of fixation
            srect = [sx1-smin sy1-smin sx1+smin sy1+smin];
            srect = srect + [-squarejitter -squarejitter squarejitter squarejitter]; % add size jitter
            rrect = [sx1+smin-(ppd*.75) sy1-rmin sx1+smin-(ppd*.25) sy1+rmin];
            rrect = rrect + [-rectjitter -rectjitter rectjitter rectjitter]; % add size jitter
        case 4
            % central fixation, large square offset, bar left of fixation
            srect = [sx1-smin sy1-smin-(sgrowth*duration) sx1+smin+(2*sgrowth*duration) sy1+smin+(sgrowth*duration)];
            srect = srect + [-squarejitter -squarejitter squarejitter squarejitter]; % add size jitter
            rrect = [sx1-smin+(ppd*.25) sy1-rmin sx1-smin+(ppd*.75) sy1+rmin];
    end
    
    srect1=srect;
    rrect1=rrect;
    
    Screen(sqtype,w,scolor,srect);                   %draw square
    if strcmp(bartype,'DrawTexture')
        Screen(bartype,w,bartexid,[],rrect);
    else
        Screen(bartype,w,rcolor,rrect);                         %draw bar
    end
    Screen('FillOval',w,fcolor,...                           %draw fixation
        [sx-fr,sy-fr,sx+fr,sy+fr]);
    
    arrow_points = [sx+fr+4 sy; sx+fr+4+20 sy; sx+fr+4+15 sy-5; sx+fr+4+20 sy; sx+fr+4+15 sy+5; sx+fr+4+20 sy]; % rightward arrow
    switch condition                                            %pursuit draw cue
        case {1 3 4 10}
            % right
            Screen('DrawLine',w,fcolor,arrow_points(1,1),arrow_points(1,2),arrow_points(2,1),arrow_points(2,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(3,1),arrow_points(3,2),arrow_points(4,1),arrow_points(4,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(5,1),arrow_points(5,2),arrow_points(6,1),arrow_points(6,2),3);
        case {11 13}
            % left
            arrow_points = RotatePointsPTB(arrow_points,180,[sx sy]);
            Screen('DrawLine',w,fcolor,arrow_points(1,1),arrow_points(1,2),arrow_points(2,1),arrow_points(2,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(3,1),arrow_points(3,2),arrow_points(4,1),arrow_points(4,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(5,1),arrow_points(5,2),arrow_points(6,1),arrow_points(6,2),3);
        case {2 12}
            % circle (no pursuit)
            Screen('FrameOval',w,fcolor,[sx-fr-4,sy-fr-4,sx+fr+4,sy+fr+4],2);
        case 9
            % 45 deg down/left
            arrow_points = RotatePointsPTB(arrow_points,-135,[sx sy]);
            Screen('DrawLine',w,fcolor,arrow_points(1,1),arrow_points(1,2),arrow_points(2,1),arrow_points(2,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(3,1),arrow_points(3,2),arrow_points(4,1),arrow_points(4,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(5,1),arrow_points(5,2),arrow_points(6,1),arrow_points(6,2),3);
        case -9
            % 45 deg down/right
            arrow_points = RotatePointsPTB(arrow_points,-45,[sx sy]);
            Screen('DrawLine',w,fcolor,arrow_points(1,1),arrow_points(1,2),arrow_points(2,1),arrow_points(2,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(3,1),arrow_points(3,2),arrow_points(4,1),arrow_points(4,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(5,1),arrow_points(5,2),arrow_points(6,1),arrow_points(6,2),3);            
        case {0 22.5 45 67.5 90}
            % X deg down/right (0 right, 90 down)
            arrow_points = RotatePointsPTB(arrow_points,-condition,[sx sy]);
            Screen('DrawLine',w,fcolor,arrow_points(1,1),arrow_points(1,2),arrow_points(2,1),arrow_points(2,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(3,1),arrow_points(3,2),arrow_points(4,1),arrow_points(4,2),3);
            %Screen('DrawLine',w,fcolor,arrow_points(5,1),arrow_points(5,2),arrow_points(6,1),arrow_points(6,2),3);              
    end
    
    lastflip = Screen('Flip',w);
    
    WaitSecs(fixpause);
    
    miss=0;
    if eyetracking==1
        result=iViewXComm('send', ivx, 'ET_STR');           %datastreamingon
    end
    
    for i=0:duration
        
        
        switch condition                                    %animation parameters
            case 1
                % fix pursuit, sq grow+trans, rect static (+growth)
                srect(2)=srect1(2)-(sgrowth*i);
                srect(3)=srect1(3)+(2*sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i);
                sx=sx1+(sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i);
            case 2
                % fix static, square grows, rect trans_away_fix (+growth)
                srect(1)=srect1(1)-(sgrowth*i);
                srect(2)=srect1(2)-(sgrowth*i);
                srect(3)=srect1(3)+(sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i);
                rrect(1)=rrect1(1)-(sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i);
                rrect(3)=rrect1(3)-(sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i);
            case 3
                % fix pursuit, sq grow+trans, rect trans (+growth)
                sx=sx1+(sgrowth*i);
                srect(2)=srect1(2)-(sgrowth*i);
                srect(3)=srect1(3)+(2*sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i);
                rrect(1)=rrect1(1)+(sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i);
                rrect(3)=rrect1(3)+(sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i);
            case 4
                % fix pursuit, sq static, rect static (+growth)
                sx=sx1+(sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i);
            case 9
                % fix pursuit (down-left), sq grow, rect follows sq edge (+growth)
                sx=sx1-(sgrowth*i);
                sy=sy1+(sgrowth*i);
                srect(1)=srect1(1)-(2*sgrowth*i);
                srect(4)=srect1(4)+(2*sgrowth*i);
                rrect(1)=rrect1(1)-(2*sgrowth*i);
                rrect(3)=rrect1(3)-(2*sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i)+(sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i)+(sgrowth*i);
            case -9
                % fix pursuit (down-right), sq grow, rect follows sq edge (+growth)
                sx=sx1+(sgrowth*i);
                sy=sy1+(sgrowth*i);
                srect(3)=srect1(3)+(2*sgrowth*i);
                srect(4)=srect1(4)+(2*sgrowth*i);
                rrect(1)=rrect1(1)+(2*sgrowth*i);
                rrect(3)=rrect1(3)+(2*sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i)+(sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i)+(sgrowth*i);
            case 10
                % fix pursuit toward, sq grow + trans, rect follows sq edge (+growth)
                sx=sx1+(sgrowth*i);
                %srect(1)=srect1(1)-(2*sgrowth*i);
                srect(2)=srect1(2)-(sgrowth*i);
                srect(3)=srect1(3)+(2*sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i);
                rrect(1)=rrect1(1)+(2*sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i);
                rrect(3)=rrect1(3)+(2*sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i);
            case 11
                % fix pursuit away (1/2), sq grow + trans (1/2) with pursuit, rect follows sq edge (+growth)
                sx=sx1-(.5*sgrowth*i);
                srect(1)=srect1(1)-(1.5*sgrowth*i);
                srect(2)=srect1(2)-(sgrowth*i);
                srect(3)=srect1(3)+(.5*sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i);
                rrect(1)=rrect1(1)+(.5*sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i);
                rrect(3)=rrect1(3)+(.5*sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i);
            case 12
                % fix static, square grows, rect trans_away_fix (+growth)
                srect(1)=srect1(1)-(sgrowth*i);
                srect(2)=srect1(2)-(sgrowth*i);
                srect(3)=srect1(3)+(sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i);
                rrect(1)=rrect1(1)+(sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i);
                rrect(3)=rrect1(3)+(sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i);
            case 13
                % fix pursuit away, sq grow, rect static (+growth)
                sx=sx1-(sgrowth*i);
                srect(1)=srect1(1)-(2*sgrowth*i);
                srect(2)=srect1(2)-(sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i);
            case {0 22.5 45 67.5 90}
                % fix pursuit (down-right), sq grow, rect follows sq edge (+growth)
                % dx,dy defines the pursuit vector
                switch condition
                    case 0
                        dx = 1;
                        dy = 0;
                    case 22.5
                        dx = 1;
                        dy = 0.5;
                    case 45
                        dx = 1;
                        dy = 1;
                    case 67.5
                        dx = 0.5;
                        dy = 1;
                    case 90
                        dx = 0;
                        dy = 1;
                end
                
                % orig pos + pursuit translation
                sx=sx1+(dx*sgrowth*i);
                sy=sy1+(dy*sgrowth*i);
                
                
                % orig pos + box growth + pursuit translation
                srect(1)=srect1(1)-(sgrowth*i)+(dx*sgrowth*i);
                srect(2)=srect1(2)-(sgrowth*i)+(dy*sgrowth*i);
                srect(3)=srect1(3)+(sgrowth*i)+(dx*sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i)+(dy*sgrowth*i);

                % orig pos + box growth + pursuit translation
                rrect(1)=rrect1(1)+(sgrowth*i)+(dx*sgrowth*i);
                rrect(3)=rrect1(3)+(sgrowth*i)+(dx*sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i)+(dy*sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i)+(dy*sgrowth*i);
                
            case {100 122.5 145 167.5 190}
                % fix pursuit (down-right), sq grow, rect follows sq edge (+growth)
                % dx,dy defines the pursuit vector
                ang = condition - 100; % convert condition back to angle from horizontal
                dx = 1; dy = 0; % vector for horizontal pursuit
                dxy = RotatePointsPTB([dx dy],-ang,[0 0]);
                
                % orig pos + pursuit translation
                sx=sx1+(dx*sgrowth*i);
                sy=sy1+(dy*sgrowth*i);
                
                % orig pos + box growth + pursuit translation
                srect(1)=srect1(1)-(sgrowth*i)+(dx*sgrowth*i);
                srect(2)=srect1(2)-(sgrowth*i)+(dy*sgrowth*i);
                srect(3)=srect1(3)+(sgrowth*i)+(dx*sgrowth*i);
                srect(4)=srect1(4)+(sgrowth*i)+(dy*sgrowth*i);
                
                % orig pos + box growth + pursuit translation
                rrect(1)=rrect1(1)+(sgrowth*i)+(dx*sgrowth*i);
                rrect(3)=rrect1(3)+(sgrowth*i)+(dx*sgrowth*i);
                rrect(2)=rrect1(2)-(rgrowth*i)+(dy*sgrowth*i);
                rrect(4)=rrect1(4)+(rgrowth*i)+(dy*sgrowth*i);
                
            otherwise
                Priority(0);                                    %reallocate memory
                ShowCursor;
                ListenChar(0);                                  %unsurpress keyboard output
                Screen('CloseAll');
                error('animation sequence not defined for condition %d',condition)                
        end
        
        Screen(sqtype,w,scolor,srect);                   %draw square
        if strcmp(bartype,'DrawTexture')
            Screen(bartype,w,bartexid,[],rrect);
        else
            Screen(bartype,w,rcolor,rrect);                         %draw bar
        end
        Screen('FillOval',w,fcolor,...                           %draw fixation
            [sx-fr,sy-fr,sx+fr,sy+fr]);
        
        xpoint=sx;
        ypoint=sy;
        
        % the additional code for thisflip and lastflip is used to check
        % the timing of execution.  ideally, each iteration of this loop
        % should last for 1 frame, so the time between flips should be
        % equal to 1/refresh
        thisflip = Screen('Flip',w);
        if checktiming
            fprintf('Time between flips = %.04f\n',thisflip-lastflip);
            lastflip = thisflip;
        end
        
        if eyetracking==1
            
            [xeye_new, yeye_new]=EyeTrackerGazeCheck(ivx,xeye_cur,yeye_cur,w);                      %track eye position
            xeye_cur = xeye_new;
            yeye_cur = yeye_new;
            xdiff = xpoint-xeye_new;
            ydiff = ypoint-yeye_new;
            delta=sqrt((xdiff^2)+(ydiff^2));
            
            %Add this to discount trials where fixation is broken
            %             if delta>acceptable_error
            %                 miss=miss+1;
            %                 if condition ==2
            %                     if miss>acceptable_nofix 
            %                     
            %                         botch=1;
            %                     end
            %                 end
            %             end
            
            Pos_v_Trace(z,1,i+1)=xpoint;                                  %record eye position versus target position on each screen refresh
            Pos_v_Trace(z,2,i+1)=ypoint;
            Pos_v_Trace(z,3,i+1)=xeye_new;
            Pos_v_Trace(z,4,i+1)=yeye_new;
            Pos_v_Trace(z,5,i+1)=xdiff;
            Pos_v_Trace(z,6,i+1)=ydiff;
            Pos_v_Trace(z,7,i+1)=delta;
            Pos_v_Trace(z,8,i+1)=miss;
        else
            Pos_v_Trace=0;
        end
    end
    
    if eyetracking==1
        result=iViewXComm('send', ivx, 'ET_EST'); %datastreamingoff
    end
    
    [keyisdown, secs, keycode] = KbCheck(dev_ID);
    while 1                                                     %leave fixation spot at last location till answer is given
        [keyisdown, secs, keycode] = KbCheck(dev_ID);
        
        Screen('FillOval',w,fcolor,...                           %draw fixation
            [sx-fr,sy-fr,sx+fr,sy+fr]);
        
        Screen('Flip',w);
        
        if keycode(a)                   %record answer
            rawdata(z,5)=1;
            break
        elseif keycode(k)
            rawdata(z,5)=2;
            break
        end    
    end
    
    if botch==0;
        rawdata(z,8)=1;                 %mark as succesful trial
        trial_list(thistrial)=0;           %remove trial from possible trials list
        com_list=[com_list;thistrial];  %add to com_list
        
        if sc_list(cidx,sidx,iidx)<=length(growth_list)/2       %if shrinking
            if rawdata(z,5)==1;                                 %if answer shrinking
                sc_list(cidx,sidx,iidx)=sc_list(cidx,sidx,iidx)+1;  %move it in growing direction
                if sc_list(cidx,sidx,iidx)>length(growth_list);     %don't go past maximum growth value
                    sc_list(cidx,sidx,iidx)=length(growth_list);
                end
                if prevan_list(cidx,sidx,iidx)==0;                  %if previous answer was incorrect record reversal
                    reversal_list(cidx,sidx,iidx,place_list(cidx,sidx,iidx))=growth_list(sc_list(cidx,sidx,iidx));
                    place_list(cidx,sidx,iidx)=place_list(cidx,sidx,iidx)+1;
                    rawdata(z,7)=1;
                else
                    rawdata(z,7)=0;
                end
                rawdata(z,6)=1;             %mark as correct
            else                            %if anser is growing
                sc_list(cidx,sidx,iidx)=sc_list(cidx,sidx,iidx)-1;   %move in shrinking direction
                if sc_list(cidx,sidx,iidx)<1                            %don't go past minimum growth value
                    sc_list(cidx,sidx,iidx)=1;
                end
                if prevan_list(cidx,sidx,iidx)==1;                      %if previous answer was correct record reversal
                    reversal_list(cidx,sidx,iidx,place_list(cidx,sidx,iidx))=growth_list(sc_list(cidx,sidx,iidx));
                    place_list(cidx,sidx,iidx)=place_list(cidx,sidx,iidx)+1;
                    rawdata(z,7)=1;
                else
                    rawdata(z,7)=0;
                end
                rawdata(z,6)=0;         %mark as incorrect
            end
        else   %if growing
            if rawdata(z,5)==2;  %if answer is growing
                sc_list(cidx,sidx,iidx)=sc_list(cidx,sidx,iidx)-1;      %move it in shrinking direction
                if sc_list(cidx,sidx,iidx)<1                            %min value
                    sc_list(cidx,sidx,iidx)=1;
                end
                if prevan_list(cidx,sidx,iidx)==0;                      %reversal
                    reversal_list(cidx,sidx,iidx,place_list(cidx,sidx,iidx))=growth_list(sc_list(cidx,sidx,iidx));
                    place_list(cidx,sidx,iidx)=place_list(cidx,sidx,iidx)+1;
                    rawdata(z,7)=1;
                else
                    rawdata(z,7)=0;
                end
                rawdata(z,6)=1;
            else                %if answer is shrinking
                sc_list(cidx,sidx,iidx)=sc_list(cidx,sidx,iidx)+1;  %move in growing direction
                if sc_list(cidx,sidx,iidx)>length(growth_list);     %max value
                    sc_list(cidx,sidx,iidx)=length(growth_list);
                end
                if prevan_list(cidx,sidx,iidx)==1;                  %reversal
                    reversal_list(cidx,sidx,iidx,place_list(cidx,sidx,iidx))=growth_list(sc_list(cidx,sidx,iidx));
                    place_list(cidx,sidx,iidx)=place_list(cidx,sidx,iidx)+1;
                    rawdata(z,7)=1;
                else
                    rawdata(z,7)=0;
                end
                rawdata(z,6)=0;
            end
        end
    else
        rawdata(z,8)=0;  %mark as unsuccesful trial for broken fixation
    end
    save(datafile,'rawdata','reversal_list','Pos_v_Trace');      %save trial and response data
    fixcount=0;
    while KbCheck(dev_ID);end
    
end

while KbCheck(dev_ID); end

if eyetracking==1
    [result, ivx]=iViewXComm('close', ivx);   %close the connection
end

Screen('TextSize',w,24);                            %Exit text
text='Thank you for your participation!';
width=RectWidth(Screen('TextBounds',w,text));
Screen('DrawText',w,text,xc-width/2,yc,tcolor);
text='Please let the experimenter know you are done.';
width=RectWidth(Screen('TextBounds',w,text));
Screen('DrawText',w,text,xc-width/2,yc+50,tcolor);
Screen('Flip',w);

KbWait(dev_ID);

save(datafile_full);                            %save trial and response data
save(datafile,'rawdata','reversal_list','Pos_v_Trace');
Priority(0);                                    %reallocate memory
ShowCursor;
ListenChar(0);                                  %unsurpress keyboard output
Screen('CloseAll');
