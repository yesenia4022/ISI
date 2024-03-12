function makeRainTexture_beta

global Mstate screenPTR screenNum 

global GtxtrAll OriAll StimLoc Dim_source Dim_dest  %'playgrating' will use these

Screen('Close')  %First clean up: Get rid of all textures/offscreen windows

GtxtrAll = []; OriAll = []; StimLoc = []; Dim_dest = []; Dim_source = [];  %reset

screenRes = Screen('Resolution',screenNum);
%Replace it with this (comes from flipInterval)
screenRes.hz = Mstate.refresh_rate;

pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

P = getParamStruct;

barW = 2*pi*Mstate.screenDist*P.barWidth/360;  %bar width in cm
barL = 2*pi*Mstate.screenDist*P.barLength/360;  %bar length in cm


%Make black/white domain for texture creation (these are different from the bw/color domains created below)
if P.bw_bit == 0
    bwdom = -1;
elseif P.bw_bit == 1
    bwdom = 1;
else
    bwdom = [-1 1];
end

colordom = getColorDomain(P.colorspace);
%%%%%%%%%%%%%%%%%%

%Make orientation domain (also used later)
if P.speed ~= 0
    oridom = linspace(P.ori,P.ori+360,P.n_ori+1);  %It goes to 360 because it is actually 'direction'
else 
    oridom = linspace(P.ori,P.ori+180,P.n_ori+1);
end
oridom = oridom(1:end-1);


Dim_dest = [round(pixpercmY*barL) round(pixpercmX*barL)]; %Square on the screen
Dim_source = [Dim_dest(1) Dim_dest(1)]/P.x_zoom; %Size of texture written to card. Keep square for ease.

Im = zeros(Dim_source); 
barWpix = round(barW*pixpercmY/P.x_zoom);  %width in pixels of source texture (not the destination)
id1 = round(size(Im,2)/2 - barWpix/2)+1;
id2 = id1+barWpix;
Im(:,id1:id2) = 1;
Im = Im*P.contrast/100;

for oriid = 1:length(oridom)
    Imdum = imrotate(Im,oridom(oriid),'nearest','crop');
    for bwid = 1:length(bwdom)
        for colorid = 1:length(colordom)
            Gtxtr(oriid,bwid,colorid) = putinTexture(Imdum*bwdom(bwid),colordom(colorid),P);
        end
    end
end

%%%
%%This next part used to be in 'playrain', but for really big stimuli it took too long and messed up COM timing%%%%%
%%%

%The following assumes the screen is curved
xcm = 2*pi*Mstate.screenDist*P.x_size/360;  %stimulus width in cm
xN = round(xcm*pixpercmX);  %stimulus width in pixels
ycm = 2*pi*Mstate.screenDist*P.y_size/360;   %stimulus height in cm
yN = round(ycm*pixpercmY);  %stimulus height in pixels

% xcm = 2*Mstate.screenDist*tan(P.x_size*pi/180/2);  %stimulus width in cm
% xN = round(xcm*pixpercmX)  %stimulus width in pixels
% ycm = 2*Mstate.screenDist*tan(P.y_size*pi/180/2);   %stimulus height in cm
% yN = round(ycm*pixpercmY);  %stimulus height in pixels

%These define the perimeters of the "location grid"
xran = [P.x_pos-ceil(xN/2)+1  P.x_pos+floor(xN/2)];
yran = [P.y_pos-ceil(yN/2)+1  P.y_pos+floor(yN/2)];

%%%Make xy domains%%%
xdom = linspace(xran(1),xran(2),P.Nx); %Make x domain  (these are the center locations of the bar)
ydom = linspace(yran(1),yran(2),P.Ny); %Make y domain

%Make bw domain
if P.bw_bit == 0 || P.bw_bit == 1
    bwdom = 1;
else
    bwdom = [1 2];
end


N_Im = round(P.stim_time*screenRes.hz/P.h_per); %number of images to present

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create independent sequence for each parameter

%Preallocate!!!
StimLoc = cell(N_Im,P.h_per);
OriAll = cell(1,N_Im);
GtxtrAll = cell(1,N_Im);
for i = 1:N_Im
    for j = 1:P.h_per
        StimLoc{i,j} = zeros(4,P.Ndrops);        
    end
    OriAll{i} = zeros(1,P.Ndrops);
    GtxtrAll{i} = zeros(1,P.Ndrops);
end

[xmesh ymesh orimesh bwmesh colormesh] = ndgrid(1:length(xdom),1:length(ydom),1:length(oridom),1:length(bwdom),1:length(colordom));

for k = 1:P.Ndrops
    
    seqdum = [];
    for l = 1:ceil(N_Im/length(xmesh(:)))
        %Used this one until 8/10/18.  It is more periodic and leaves holes.
        %s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed+(l-1)*100 + (k-1)*1000);
        s = RandStream('mt19937ar','Seed',P.rseed*100 + (l-1)*100 + (k-1)*1000);
        [dumval dumid] = sort(rand(s,1,length(xmesh(:))));
        seqdum = [seqdum dumid];
    end
    
    xseqID(:,k) = xmesh(seqdum(1:N_Im));  %The sequence for each drop
    yseqID(:,k) = ymesh(seqdum(1:N_Im));  %The sequence for each drop
    oriseqID(:,k) = orimesh(seqdum(1:N_Im));  %I don't know why, but this one should be an index, not the actual ori
    bwseqID(:,k) = bwmesh(seqdum(1:N_Im));  %The sequence for each drop
    colorseqID(:,k) = colormesh(seqdum(1:N_Im));  %The sequence for each drop
    
    %How I used to do it (pre 3/23/2016)
%     xseqId = round(rand(s,1,N_Im)*length(xdom)+.5);  %The sequence for each drop
%     xseq = xdom(xseqId);
% 
%     yseqId = round(rand(s,1,N_Im)*length(ydom)+.5);
%     yseq = ydom(yseqId);
% 
%     oriseq = round(rand(s,1,N_Im)*length(oridom)+.5);

 
%     figure,hist(xseq,length(xdom))
%     figure,hist(yseq)
%     figure,hist(oriseq,length(oridom))
%     figure,hist(colorseq,length(colordom))
%     figure,hist(bwseq,length(bwdom))

    oriseq = oridom(oriseqID(:,k))';
    xseq = xdom(xseqID(:,k))';
    yseq = ydom(yseqID(:,k));
    
    oriseq = oriseq(:);
    xseq = xseq(:);
    yseq = yseq(:);

    if strcmp(P.gridType,'polar')

        xseqRot = (xseq-P.x_pos).*cos(-oriseq*pi/180) - (yseq-P.y_pos).*sin(-oriseq*pi/180);
        yseqRot = (xseq-P.x_pos).*sin(-oriseq*pi/180) + (yseq-P.y_pos).*cos(-oriseq*pi/180);

        yseqRot = round(yseqRot*pixpercmY/pixpercmX); %correct for non-square pixels
        
        xseq = xseqRot+P.x_pos; %Center location of texture on the screen, in pixel units
        yseq = yseqRot+P.y_pos;
 
    end

%     bwseq = round(rand(s,1,N_Im)*length(bwdom)+.5); %this should remain an index value
%     
%     colorseq = round(rand(s,1,N_Im)*length(colordom)+.5); %this should remain an index value

    xseqL = xseq-(ceil(Dim_dest(2)/2)-1); %Destination window, in pixels
    xseqR = xseq+floor(Dim_dest(2)/2);
    yseqL = yseq-(ceil(Dim_dest(1)/2)-1);
    yseqR = yseq+floor(Dim_dest(1)/2);

    Dinc = 2*Mstate.screenDist*tan(P.speed/2*pi/180);  %cm increment per frame

    for i = 1:N_Im
        xinc = Dinc*cos(oriseq(i)*pi/180);
        yinc = -Dinc*sin(oriseq(i)*pi/180);  %negative because origin is at top
        for j = 1:P.h_per
            dx = (j-1)*xinc;
            dx = round(dx*pixpercmX);  %convert to pixels
            dy = (j-1)*yinc;
            dy = round(dy*pixpercmY);  %convert to pixels
            xseqL2 = xseqL(i)+dx;
            xseqR2 = xseqR(i)+dx;
            yseqL2 = yseqL(i)+dy;
            yseqR2 = yseqR(i)+dy;

            StimLoc{i,j}(:,k) = [xseqL2 yseqL2 xseqR2 yseqR2]';
        end

        GtxtrAll{i}(k) = Gtxtr(oriseqID(i,k),bwseqID(i,k),colorseqID(i,k));


    end

end

domains = struct;
domains.oridom = oridom;
domains.xdom = xdom;
domains.ydom = ydom;
domains.bwdom = bwdom;
domains.colordom = colordom;

if Mstate.running %if its in the looper
    
    Pseq.oriseq = oriseqID;
    Pseq.xseq = xseqID;
    Pseq.yseq = yseqID;
    Pseq.bwseq = bwseqID;
    Pseq.colorseq = colorseqID;
    
    saveLog_rain(domains,Pseq)  %append log file with the latest sequence

end




function Gtxtr = putinTexture(Im,colortype,P)

global screenPTR

%%%%%%%%%%%%%%%%%%%%%%%
%Equate Contrast: This is a total hack%%
if strcmp(P.colorspace,'DKL')
    switch colortype
        case 4 %S
            Im = Im*.15/.82 * 3;
        case 5 %L-M
            Im = Im;
        case 6 %L+M
            Im = Im*.15/1.0;
    end
elseif strcmp(P.colorspace,'LMS')
    switch colortype
        case 2 %L
            Im = Im;
        case 3 %M
            Im = Im*.2/.23;
        case 4 %S
            Im = Im*.2/.82 * 3;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imcell{1} = Im;
Idraw = ImtoRGB(Imcell,colortype,P,[]);
Gtxtr = Screen(screenPTR, 'MakeTexture', Idraw);