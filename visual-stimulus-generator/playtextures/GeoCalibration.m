function GeoCalibration

%Before running:
% Set the screen from the Master Main Window, but make sure that updateMon at the
% Display is not using a geocalibration file. Save calibration file to a
% location pointed to by updateMon flip uv projectors orientation

global Mstate screenPTR screenNum 

global Gtxtr TDim xyin xyout   %Created in makeGratingTexture

global imcolortrx   %This should be unity during calibration

P = getParamStruct;

%Pstruct = getParamStruct;

SquareW = P.DotSize; %Size in pixels on screen
DotSpacing = P.Dot_Spacing; %pixels

screenRes = Screen('Resolution',screenNum);
pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

Screen(screenPTR, 'FillRect', 0); %Set Background

HideCursor

if ~P.calib_check
    
    xyin = [];
    xyout = [];
  
    %initialize the texture
    dum = zeros(1,1,3); dum(2) = 128; %dum(1) = 255;
    Gtxtr = Screen(screenPTR, 'MakeTexture',dum);
    dum = zeros(1,1,3); dum(2) = 255;
    Gtxtr(2) = Screen(screenPTR, 'MakeTexture',dum); %He will be highlighted to indicate that its his turn
    dum = zeros(1,1,3); dum(3) = 255; % Blue gun is the UV. 
    Gtxtr(3) = Screen(screenPTR, 'MakeTexture',dum);  %This will be at the tip of the mouse pointer
    
    TDim = [SquareW SquareW];
    StimPiece = [0 0 TDim(2)-1 TDim(1)-1]';
    
    marg = round(DotSpacing/2);
    marg = 100;
    xdim = (marg:DotSpacing:(screenRes.width-marg))+25;
    ydim = (marg:DotSpacing:(screenRes.height-marg))+100;
    [x y] = meshgrid(xdim,ydim);
    
    Ndots = length(x(:));
    
    %STCentered = StimPiece - [TDim(2)/2 TDim(1)/2 TDim(2)/2 TDim(1)/2];
    StimLocs = StimPiece*ones(1,length(x(:)));
    StimLocs([1 3],:) = StimLocs([1 3],:)+[x(:)'; x(:)'];
    StimLocs([2 4],:) = StimLocs([2 4],:)+[y(:)'; y(:)'];
    
    StimPieces = StimPiece*ones(1,length(x(:)));
    %
    % Screen('DrawTextures', screenPTR,Gtxtr(1),StimPieces(:,2:end),StimLocs(:,2:end));
    % Screen('DrawTextures', screenPTR,Gtxtr(2),StimPieces(:,1),StimLocs(:,1));
    % Screen('Flip', screenPTR);
    
    bLast = [0 0 0];
    keyIsDown = 0;
    i = 1;
    while ~keyIsDown & i<=Ndots
        
        [mx,my,b] = GetMouse(screenPTR);

        db = bLast - b; %'1' is a button release
        
        %%%Left Button release%%%
        if ~sum(abs([1 0 0]-db))
            
            mx_all(i) = mx;
            my_all(i) = my;
            
            i = i+1;
            
        end
        
        if i<=Ndots %if its not about to be the last one
            lowdom = 1:Ndots;
            lowdom(i) = []; %index of next dot that should be highlighted

            Screen('DrawTextures', screenPTR,Gtxtr(1),StimPieces(:,lowdom),StimLocs(:,lowdom));
            Screen('DrawTextures', screenPTR,Gtxtr(2),StimPieces(:,i),StimLocs(:,i));
            
            mouseSquareLoc = [0 0 TDim(2)-1 TDim(1)-1]' + [mx my mx my]' ;
            Screen('DrawTextures', screenPTR,Gtxtr(3),StimPieces(:,1),mouseSquareLoc);
            
            Screen('Flip', screenPTR);
        end
        
        bLast = b;
        
        keyIsDown = KbCheck();
        
    end
    
    xyin = [x(:) y(:)];
    xyout = [mx_all' my_all'];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Screen(screenPTR, 'FillRect', 128)
    Screen(screenPTR, 'Flip');
    
    Screen('Close')  %Get rid of all textures/offscreen windows
    sca
    uisave({'xyin','xyout'},'xyinout')
    %save xyinout xyin xyout
    
else
    
    ch = 1; %channel to undergo calib transformation
    
    %initialize the texture
    dum = zeros(1,1,3); dum(2) = 128; %green dots
    Gtxtr = Screen(screenPTR, 'MakeTexture',dum);
    dum = zeros(1,1,3); dum(ch) = 255; %UV dots
    Gtxtr(2) = Screen(screenPTR, 'MakeTexture',dum);
    
    TDim = [SquareW SquareW];
    StimPiece = [0 0 TDim(2)-1 TDim(1)-1]';
    
    xdim = DotSpacing/2:DotSpacing:screenRes.width;
    ydim = DotSpacing/1.5:DotSpacing:screenRes.height;
    [x y] = meshgrid(xdim,ydim);
        
    oneMat = ones(size(y));
    xUV = imcolortrx{ch}(1,1)*x + imcolortrx{ch}(2,1)*y + imcolortrx{ch}(3,1)*y.*x + imcolortrx{ch}(4,1)*oneMat;
    yUV = imcolortrx{ch}(1,2)*x + imcolortrx{ch}(2,2)*y + imcolortrx{ch}(3,2)*y.*x + imcolortrx{ch}(4,2)*oneMat;
    
    Ndots = length(x(:));
    
    StimLocs = StimPiece*ones(1,length(x(:)));
    StimLocs(1,:) = StimLocs(1,:)+x(:)';
    StimLocs(3,:) = StimLocs(3,:)+x(:)';
    StimLocs(2,:) = StimLocs(2,:)+y(:)';
    StimLocs(4,:) = StimLocs(4,:)+y(:)';
    
    StimLocsUV = StimPiece*ones(1,length(x(:)));
    StimLocsUV(1,:) = StimLocsUV(1,:)+xUV(:)';
    StimLocsUV(3,:) = StimLocsUV(3,:)+xUV(:)';
    StimLocsUV(2,:) = StimLocsUV(2,:)+yUV(:)';
    StimLocsUV(4,:) = StimLocsUV(4,:)+yUV(:)';
    
    StimPieces = StimPiece*ones(1,length(x(:)));
    
    
    keyIsDown = 0; 
    while ~keyIsDown        
        
        Screen('DrawTextures', screenPTR,Gtxtr(1),StimPieces,StimLocs);
        Screen('DrawTextures', screenPTR,Gtxtr(2),StimPieces,StimLocsUV);
        Screen('Flip', screenPTR);
        
        keyIsDown = KbCheck();
        
    end
    
    Screen('Close')  %Get rid of all textures/offscreen windows
    sca
    
end

