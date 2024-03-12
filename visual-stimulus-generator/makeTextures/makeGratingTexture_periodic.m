function makeGratingTexture_periodic

%make one cycle of the grating

global Mstate screenPTR screenNum %movieBlock 

global Gtxtr TDim  %'playgrating' will use these

global imcolortrx %set in updateMonitor.m

Screen('Close')  %First clean up: Get rid of all textures/offscreen windows

Gtxtr = []; TDim = [];  %reset

% frame1=1;

[P] = getParamStruct;
screenRes = Screen('Resolution',screenNum);

pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;


if strcmp(P.altazimuth,'none')
    
    %The following assumes the screen is curved
    xcm = 2*pi*Mstate.screenDist*P.x_size/360;  %stimulus width in cm
    xN = round(xcm*pixpercmX);  %stimulus width in pixels
    ycm = 2*pi*Mstate.screenDist*P.y_size/360;   %stimulus height in cm
    yN = round(ycm*pixpercmY);  %stimulus height in pixels
    
else
    
    %The following assumes a projection of spherical coordinates onto the
    %flat screen
    xN = 2*Mstate.screenDist*tan(P.x_size/2*pi/180);  %grating width in cm
    xN = round(xN*pixpercmX);  %grating width in pixels
    yN = 2*Mstate.screenDist*tan(P.y_size/2*pi/180);  %grating height in cm
    yN = round(yN*pixpercmY);  %grating height in pixels
    
end

xN = round(xN/P.x_zoom);  %Downsample for the zoom
yN = round(yN/P.y_zoom);
if ~yN
    yN = 1;% For the case of blanks
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maskAZM = P.altazimuth2;
if P.ori2 == 0 | P.ori2 == 180;
    maskAZM = 'azimuth';
elseif P.ori2 == 90 | P.ori2 == 270
    maskAZM = 'altitude';
end

carrierAZM = P.altazimuth;
if P.ori == 0 | P.ori == 180
    carrierAZM = 'azimuth';
elseif P.ori == 90 | P.ori == 270
    carrierAZM = 'altitude';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%

%create the mask
xdom = linspace(-P.x_size/2,P.x_size/2,xN);
ydom = linspace(-P.y_size/2,P.y_size/2,yN);
[xdom ydom] = meshgrid(xdom,ydom);
r = sqrt(xdom.^2 + ydom.^2);
if strcmp(P.mask_type,'disc')
    mask = zeros(size(r));
    id = find(r<=P.mask_radius);
    mask(id) = 1;
elseif strcmp(P.mask_type,'gauss')
    mask = exp((-r.^2)/(2*P.mask_radius^2));
elseif strcmp(P.mask_type,'grating')
    %this mask is defined by the plaid parameters
    
%     [sdom_mask dum x_ecc_mask y_ecc_mask] = makeGraterDomain_beta(xN,yN,P.ori2,P.s_freq2,1,maskAZM,eye(4,2));
%     [dum mask] = makeSeparableProfiles(0,sdom_mask,x_ecc_mask,y_ecc_mask,2);
%     
%     if P.mask_floor == 0
%         %mask = mask; %[-1 1]
%     elseif P.mask_floor == 0.5
%         mask = (mask+1)/2; %[0 1] this makes it a contrast modulator
%     end
    mask = [];
else
    mask = [];
end
mask = single(mask);
%%%%%%%%%

%%%%%%
%%%%%%BETA VERSION


for q = 1:length(imcolortrx) %Spatial domain for each color
    %[sdom{q} tdom x_ecc{q} y_ecc{q}] = makeGraterDomain_beta(xN,yN,P.ori,P.s_freq,P.t_period,P.altazimuth,imcolortrx{q});%orig
    [sdom{q} tdom x_ecc{q} y_ecc{q}] = makeGraterDomain_beta(xN,yN,P.ori,P.s_freq,P.t_period,P.altazimuth,imcolortrx{q});
    
    
    if P.plaid_bit
        %I am ignoring t_period2 for now, and just setting it to t_period
        %     if strcmp(P.altazimuth,'altitude')
        %         AZ2 = 'azimuth'
        %     elseif strcmp(P.altazimuth,'azimuth')
        %         AZ2 = 'altitude';
        %     end
        %[sdom2{q} tdom2 x_ecc2{q} y_ecc2{q}] = makeGraterDomain(xN,yN,P.ori2,P.s_freq2,P.t_period2,AZ2,imcolortrx{q});
        [sdom2{q} tdom2 x_ecc2{q} y_ecc2{q}] = makeGraterDomain_beta(xN,yN,P.ori2,P.s_freq2,P.t_period2,P.altazimuth2,imcolortrx{q});
        
        if length(tdom)>length(tdom2)
            repN = ceil(length(tdom)/length(tdom2));
            tdom2 = repmat(tdom2,[1 repN]);
            tdom2 = tdom2(1:length(tdom));
        elseif length(tdom2)>length(tdom)
            repN = ceil(length(tdom2)/length(tdom));
            tdom = repmat(tdom,[1 repN]);
            tdom = tdom(1:length(tdom));
        end
        
    end
    
end

flipbit = 0;

harmdom = eval(['[' P.harmonics ']']);

if P.phaseShuff
    phaseshift = 2*pi*rand(1,length(harmdom))*100;
else
    phaseshift = 0*rand(1,length(harmdom));
end
phaseshift = phaseshift./harmdom;

Ampmax = 0;
for harmID = 1:length(harmdom)
    Ampmax = Ampmax+1/harmdom(harmID);
end
Ampmax
if ~P.separable
    
    for i = 1:length(tdom)
        
        for q = 1:length(imcolortrx)
            Im = 0;
            
             for harmID = 1:length(harmdom)

                if P.speed_bit

                    Im = Im+makePerGratFrame_insep(sdom{q}+phaseshift(harmID),tdom,i,harmdom(harmID),1) / harmdom(harmID); 

                else

                    Im = Im+makePerGratFrame_insep((sdom{q}+phaseshift(harmID))*harmdom(harmID),tdom,i,1,1) / harmdom(harmID); 

                end

                    

            end        
            Im = Im/Ampmax;
            %Im= Im/max(Im(:)); %hack for ronan
            
            %figure,imagesc(Im)
            if P.plaid_bit
                Im2 = makePerGratFrame_insep(sdom2{q},tdom2,i,1,2); %[-1 1]*contrast2
                Im2 = Im2 + P.midpoint2;  %midpoint varies between -1 and 1;
                if strcmp(P.plaid_function(1:3),'mul')
                    Im = Im.*Im2;
                elseif strcmp(P.plaid_function(1:3),'add')
                    Im = Im+Im2;
                end
            end
            
            
            if P.noise_bit
                if rem(i,P.noise_lifetime) == 1
                    %                 nwx = round(P.noise_width/P.x_zoom);
                    %                 nwy = round(P.noise_width/P.y_zoom);
                    %                 noiseIm = makeNoiseIm(size(Im),nwx,nwy,P.noise_type);
                    
                    noiseIm = makeNoiseIm_beta(size(Im),P,x_ecc{q},y_ecc{q});
                    
                    flipbit = 1-flipbit;
                    if flipbit
                        noiseIm = 1-noiseIm;
                    end
                end
                
                Im = Im - 2*noiseIm;
                Im(find(Im(:)<-1)) = -1;
                
            end
            
            Imcell{q} = Im;
            
            
        end
        
        ImRGB = ImtoRGB(Imcell,P.colormod,P,mask);
        
        Gtxtr(i) = Screen(screenPTR, 'MakeTexture', ImRGB);
        
    end
    
else
    
    for q = 1:length(imcolortrx)
        [amp{q} temp{q}] = makeSeparableProfiles(tdom,sdom{q},x_ecc{q},y_ecc{q},1);
        if P.plaid_bit
            [amp2{q} temp2{q}] = makeSeparableProfiles(tdom2,sdom2{q},x_ecc2{q},y_ecc2{q},2);
        end
    end
    
    for i = 1:length(tdom)
        
        for q = 1:length(imcolortrx)
            
            Im = amp{q}(i)*temp{q};
            
            if P.plaid_bit
                Im2 = amp2{q}(i)*temp2{q};               
                Im2 = Im2 + P.midpoint2;  %midpoint varies between -1 and 1;
                if strcmp(P.plaid_function(1:3),'mul')
                    Im = Im.*Im2;
                elseif strcmp(P.plaid_function(1:3),'add')
                    Im = Im+Im2;
                end
            end
            
            if P.noise_bit
                if rem(i,P.noise_lifetime) == 1
                    %                 nwx = round(P.noise_width/P.x_zoom);
                    %                 nwy = round(P.noise_width/P.y_zoom);
                    %                 noiseIm = makeNoiseIm(size(Im),nwx,nwy,P.noise_type);
                    
                    noiseIm = makeNoiseIm_beta(size(Im),P,x_ecc{q},y_ecc{q});
                    
                    flipbit = 1-flipbit;
                    if flipbit
                        noiseIm = 1-noiseIm;
                    end
                end
                
                %Im = Im - 2*noiseIm;
                Im = Im.*noiseIm;
                Im(find(Im(:)<-1)) = -1;
                
            end
            
            Imcell{q} = Im;
            
        end
        
        
        ImRGB = ImtoRGB(Imcell,P.colormod,P,mask);
        Gtxtr(i) = Screen(screenPTR, 'MakeTexture', ImRGB);
        
    end
    
end





TDim = size(ImRGB(:,:,1));
TDim(3) = length(Gtxtr);


