function grating_stimuli_contrast_levels(t, repetitions, dx, withgray_yn, random_yn, freqCyclesPerPix, outputfile)
% grating stimuli presented at a set of angles (defined below) with speed dx

global win texid;

%% Init variables
if ~exist('outputfile', 'var')
    outputfile = 'C:/Data/EPhys/stimuli_info/';
end

if ~exist('repetitions', 'var')
    repetitions=1;
end

if ~exist('withgray_yn', 'var')
    withgray_yn=1;
end

if ~exist('random_yn', 'var')|| isempty(random_yn)
    random_yn=1;
end

if ~exist('freqCyclesPerPix', 'var')|| isempty(freqCyclesPerPix)
    % Grating frequency in cycles / pixel
    freqCyclesPerPix = 0.01;  %this will correspond to 100pixels swept in 1 sec
end

% Drift speed cycles per second
if ~exist('dx', 'var') || isempty(dx)
    dx=5;
    cyclesPerSecond = 1;
else
    cyclesPerSecond =dx;
end

%% Set all params
winRect = Screen('Rect', win);
texRect = Screen('Rect', texid(2));

white = 1;%Screen('ColorRAnge', win);

grey = white / 2;
ifi = Screen('GetFlipInterval',win);
contrast_levels=uint8([0,255; 2,127;12,32]);
nlevels=size(contrast_levels,1);

rseed=rng;

%% syncing red mask surrounding the stimuli area
%take the area around the screen and use it show red syncing signal
maxRedIntensity = 90;
maskIm=Screen('GetImage',texid(1));
maskIm=maskIm(:,:,1);
maskIm=255-maskIm;
maskIm(maskIm<255)=0;
maskImAlpha = uint8((maskIm>0)*255);
maskImAlpha=repelem(maskImAlpha,2,1);
maskIm(maskIm==255)=maxRedIntensity;
maskIm=repelem(maskIm,2,1);
zerolayer = zeros(size(maskIm),'uint8');
maskIm=cat(3,maskIm,zerolayer,zerolayer,maskImAlpha);
redTex =  Screen('MakeTexture', win, maskIm);
%%
    
wall_offset=[0 0];

% Grating size in pixels
gratingSizePix = min(winRect(3:4))-min(wall_offset);

% Frequency in Radians
freqRad = freqCyclesPerPix * 2 * pi;

% Define Half-Size of the grating image.
texsize = gratingSizePix / 2;

% This is the visible size of the grating
% visibleSize = 2 * texsize + 1;
visibleSize = min(texRect(3)-texRect(1),texRect(4)-texRect(2));

% Make a destination rectangle for our textures and center this on the
% screen
dstRect = [0 0 visibleSize visibleSize];
dstRect_tex = CenterRect(dstRect, texRect);
dstRect = CenterRect(dstRect, winRect);


% We set PTB to wait one frame before re-drawing
waitframes = 1;
% Calculate the wait duration
waitDuration = waitframes * ifi;

% Recompute pixPerCycle, this time without the ceil() operation from above.
% Otherwise we will get wrong drift speed due to rounding errors
pixPerCycle = 1 / freqCyclesPerPix;

% Translate requested speed of the grating (in cycles per second) into
% a shift value in "pixels per frame"
shiftPerFrame = cyclesPerSecond * pixPerCycle * waitDuration;

% Define our grating. Note it is only 1 pixel high. PTB will make it a full
% grating upon drawing
scale_fact = 2;
x = meshgrid(-(texsize*scale_fact):(texsize + pixPerCycle)*2, 1);
nx_2=uint32(round(numel(x)/2));

gratingpic = (grey * cos(freqRad*x) + grey)>0.5*white;
A= repmat(gratingpic,numel(x),1);
%compute average intensity of the pattern
av_intensity=16;

% create gray texture 
clearTex=Screen('MakeTexture', win, zeros(texRect(4), texRect(3),3));
Screen('FillRect', clearTex, [0 av_intensity 0], dstRect_tex); 
Screen('DrawTexture', texid(2), clearTex);
clearTex_morphed = correct_distortion(win, texid);
Screen('Close', clearTex);


%% create the set of gratings to choose from 
% angle increment by 45 degrees
% a set consists of every angle {0:15:360}
% every grating will be presented drifting and statically
% every grating will be repeated 'repetitions' number of times
% there will be periods of gray screen before every repetition 
% black periods will be encoded as [-1,-1]
if withgray_yn
    set_st=[-1,-1];
else
    set_st=[];
end
contrast_i=1:nlevels;
alphas=0:45:359;
[aa,bb] = meshgrid(contrast_i,alphas);
allconfigs=cat(2,aa',bb');
allconfigs=reshape(allconfigs,[],2);
allconfigs_n=size(allconfigs,1);
allconfigs_reps=zeros(repetitions*allconfigs_n,1);

%% start iterating thru stimulus set
tstart = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');
tstamp = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
% Sync us to the vertical retrace
vbl = Screen('Flip', win);
% nfr=round(t*repetitions*(8+8+1)*60+(repetitions*(8+8+1)+1)*60*ipi_duration+1000);
% frame_array=zeros(684,608,3,nfr,'uint8');
% fi=1;
for ri=1:repetitions
% for i=1:nst  
    if random_yn
        cur_perm=randperm(allconfigs_n);
    else
        cur_perm=1:allconfigs_n;
    end
    allconfigs_reps((ri-1)*allconfigs_n+1:ri*allconfigs_n)=cur_perm;
    
    if withgray_yn
        frameCounter=0;
        t0=GetSecs;
        while GetSecs-t0<t %gray screen for t seconds
            Screen('DrawTexture', win, clearTex_morphed,[],winRect); 
            if mod(frameCounter,5)==0
                %activate the blending mode in order to combine two textures: red
                %frames and stimulus itself        
                Screen(win,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                Screen('DrawTexture', win, redTex);      
                Screen(win,'BlendFunction',GL_ONE, GL_ZERO);
            end   
            Screen('DrawingFinished', win);    
            % Flip to the screen on the next vertical retrace
            vbl = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi, 2); 
            frameCounter=frameCounter+1;
        end
    end
    for ni=1:allconfigs_n
        alpha =   allconfigs(cur_perm(ni),2);
        ci = allconfigs(cur_perm(ni),1);
       
        %setup grating pattern        
        %adjust the contrast
        Ac=uint8(A);
        Ac(A==0)=contrast_levels(ci,1);
        Ac(A==1)=contrast_levels(ci,2);
        %rotate the pattern
        B=imrotate(Ac,-alpha);
        cntr=uint32(round(length(B)/2));
        %Make a two layer mask filled with the background colour
        C=uint8(B(cntr-nx_2+1:cntr+nx_2-1, cntr-nx_2+1:cntr+nx_2-1)*white);
        C=cat(3,C,C,C);
        C(:,:,1)=0;
        C(:,:,3)=0;
        % Make the grating mask texture
        gratingMaskTex = Screen('MakeTexture', win, C);
       
        if mod(alpha,180)==0
            sinalpha=0;
        else
            sinalpha=sin(alpha*pi/180);
        end
        if mod(alpha,90)==0 && mod(alpha,180)~=0
            cosalpha=0;
        else
            cosalpha=cos(alpha*pi/180);
        end

        vsStart=floor((length(C)-visibleSize)/2);
        srcRect=[]; 
        frameCounter=0;
        t0=GetSecs;
        te=GetSecs-t0;
        while te < t+3 %t for moving, 2sec static before, and 1sec static after
            if te < 2 || te > t+2
                speed=0;
            else
                speed=dx;
            end
            if speed==0 && isempty(srcRect)
                 srcRect = [vsStart vsStart vsStart+visibleSize vsStart+visibleSize];
            elseif speed~=0
                % Calculate the xoffset for the window thru which to sample grating
                offset= mod(frameCounter * shiftPerFrame, pixPerCycle);
                yoffset = offset*sinalpha;
                xoffset = offset*cosalpha;
                % Define our source rectangle for grating sampling
                srcRect = [vsStart+xoffset vsStart+yoffset  vsStart+xoffset+ visibleSize vsStart+yoffset+visibleSize];                
            end      
            
            % Draw grating texture with black background
            Screen('FillRect', texid(2), [0,0,0]);
            Screen('DrawTexture', texid(2), gratingMaskTex, srcRect);
            morphtexid = correct_distortion(win, texid);              
            Screen('DrawTexture', win, morphtexid,[],winRect);            
            
            % draw red surrounding area every 5th frame
            if mod(frameCounter,5)==0  && speed~=0
                %activate the blending mode in order to combine two textures: red
                %frames and stimulus itself        
                Screen(win,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                Screen('DrawTexture', win, redTex);      
                Screen(win,'BlendFunction',GL_ONE, GL_ZERO);
            end   

            Screen('DrawingFinished', win);    
            % Flip to the screen on the next vertical retrace
            vbl = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi, 2);
            Screen('Close', morphtexid);   

			if speed~=0
                frameCounter=frameCounter+1;
            end	
			
            te=GetSecs-t0;
%             %%%
%             frame_array(:,:,:,fi)=Screen('GetImage', win);
%             fi=fi+1;
        end
        Screen('Close', gratingMaskTex);   
    end
   
end
Screen('Close', redTex);   
% v = VideoWriter('allframes.avi');
% 
% for fii=1:fi-1    
%     open(v);
%     writeVideo(v,frame_array(:,:,:,fii));
% end
% close(v);

tend = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');

%%save all params to the file
%write stats to the stimuli.log file
filestimuli_info = [outputfile, 'stimuli_', mfilename(),'_',tstamp,'_log.mat'] ;

% %save stats to the file
stim_name = mfilename();
save(filestimuli_info,'stim_name','tstart','tend','t','ifi','repetitions','dx','withgray_yn','freqCyclesPerPix','contrast_levels','allconfigs','allconfigs_reps');
