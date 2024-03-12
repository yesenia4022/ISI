function ImRGB = ImtoRGB(Im,colormod,P,mask)

%This function takes the image 'Im', which ranges from -1 to 1 and converts
%it into the RGB image with color direction based on colormod.  The output
%image can be made as a texture.

%colormod as an input is used for the flash grater, when color is
%randomized. 


switch colormod
    
    case 1  %this obeys the input gain values
        
        rgain = P.redgain; 
        ggain = P.greengain;
        bgain = P.bluegain;
        
    case 2  %L-cone isolation
        
        %monk aa6
%         rgain = 1;
%         bgain = -.162;
%         ggain = -.00547;

%         %6/20/10 (2 deg fundamentals)
%         rgain = 1;
%         ggain = -.1961;
%         bgain = -.0187;

        %6/20/10 (10 deg fundamentals)
%         rgain = 1;
%         ggain = -.1782;
%         bgain = -.0228;
        
        %9/12/10 (10 deg fundamentals)
%         rgain = 1;
%         ggain = -.1772;
%         bgain = -.0143;

        %5/9/16 (10 deg fundamentals)
%         rgain = 1;
%         ggain = -.1596;
%         bgain = -.0109;
        
        %9/20/17 CRT (10 deg fundamentals)
%         rgain = 1;
%         ggain = -.1653;
%         bgain = -.0130;

        %5/16/19 Vpixx (10 deg fundamentals)
        rgain = 1;
        ggain = -.2165;
        bgain = .0142;

        %fake cones
%         rgain = 1;
%         ggain = -.1115;
%         bgain = -.0218;
                
    case 3  %M-cone isolation
        
        %monk aa6
%         rgain = -1;
%         bgain = .470;
%         ggain = -.037;

%         %6/20/10 (2 deg fundamentals)
%         rgain = -1;
%         ggain = .5620;
%         bgain = -.0357;
        
        %6/20/10 (10 deg fundamentals)
%         rgain = -1;
%         ggain = .5142;
%         bgain = -.0229;
        
        %9/12/10 (10 deg fundamentals)
%         rgain = -1;
%         ggain = .5214;
%         bgain = -.0273

        %5/9/16 (10 deg fundamentals)
%         rgain = -1;
%         ggain = .4587;
%         bgain = -.0123;

        %9/20/17 CRT (10 deg fundamentals)
        rgain = -1;
        ggain = .4698;
        bgain = -.0105;
        
       %5/16/19 VPixx (10 deg fundamentals)
        rgain = -0.6921;
        ggain = 1;
        bgain = -.1704;

        
        %fake cones
%         rgain = -1.0;
%         ggain = .5948;
%         bgain = .0002;
        
    case 4  %S-cone isolation
        
        %monk aa6
%         rgain = .158;
%         bgain = -.196;
%         ggain = 1;

% %       %6/20/10 (2 deg fundamentals)
%         rgain = .1261;
%         ggain = -.2071;
%         bgain = 1;

%       %6/20/10 (10 deg fundamentals)
%         rgain = .2052;
%         ggain = -.2768;
%         bgain = 1;
        
        %%9/12/10 (10 deg fundamentals)
%         rgain = .2050;
%         ggain = -.2738;
%         bgain = 1;

        %%5/9/16 (10 deg fundamentals)
%         rgain = .3333;
%         ggain = -.3892;
%         bgain = 1;

        %%9/20/17 CRT (10 deg fundamentals)
%         rgain = .3427;
%         ggain = -.4080;
%         bgain = 1;
        
        
        %%5/16/19 Vpixx (10 deg fundamentals)
        rgain = .0590;
        ggain = -.4147;
        bgain = 1;
        
        %fake cones
%         rgain = .3624;
%         ggain = -.2742;
%         bgain = 1.0;
        
    case 5  %L-M 

        %monk aa6
%         rgain = 1;
%         bgain = -.3578;
%         ggain = .0257;
   
%         %6/20/10 (2 deg fundamentals)
%         rgain = 1;
%         ggain = -.3613;
%         bgain = -.0059;


        %6/20/10 (10 deg fundamentals)
%         rgain = 1;
%         ggain = -.334;
%         bgain = -.0016;

%         %9/12/10 (10 deg fundamentals)
%         rgain = 1;
%         ggain = -.3365;
%         bgain = -.0049;

        %9/20/17 CRT(10 deg fundamentals)
%         rgain = 1;
%         ggain = -.3082;
%         bgain = -0.0020;

        %%5/16/19 Vpixx (10 deg fundamentals)   38% total contrast
        rgain = 1;
        ggain = -.6852;
        bgain = 0.1028;
        
        
    case 6  %L+M

        %monk aa6
%         rgain = -1.0;
%         bgain = .9104;
%         ggain = -.0984;

%         %6/20/10 (2 deg fundamentals)       
%         rgain = .6617;
%         ggain = 1;
%         bgain = -.1803;

        %6/20/10 (10 deg fundamentals)       
%         rgain = .5107;
%         ggain = 1;
%         bgain = -.1603;
        
        
%         %9/12/10 (10 deg fundamentals)       
%         rgain = .5122;
%         ggain = 1;
%         bgain = -.1390;

        %5/9/16 (10 deg fundamentals)       
%         rgain = .1143;
%         ggain = 1;
%         bgain = -.0800;
        
        %9/20/17 CRT (10 deg fundamentals)       
%         rgain = .4635;
%         ggain = 1;
%         bgain = -.0893;

        %5/16/19 Vpixx (10 deg fundamentals)    90% contrast   
        rgain = .4927;
        ggain = 1;
        bgain = -.2047;
        
    case 7 %S + (L-M)

        %monk aa6
%         rgain = 1.0;
%         bgain = -.3782;
%         ggain = .1819;

%         %6/20/10 (2 deg fundamentals
%         rgain = 1.0;
%         ggain = -.3849;
%         bgain = .1517;

        %6/20/10 (10 deg fundamentals
%         rgain = 1.0;
%         ggain = -.3849;
%         bgain = .1517;
        
        %9/12/10 (10 deg fundamentals
        rgain = 1.0;
        ggain = -.3629;
        bgain = .1338;
        

        
        
    case 8 %S - (L-M)
 
        %monk aa6
%         rgain = -1.0;
%         bgain = .3369;
%         ggain = -.1455;

%         %6/20/10 (2 deg fundamentals)
%         rgain = -1.0;
%         ggain = .2613;
%         bgain = .0986;

        %6/20/10 (10 deg fundamentals)
%         rgain = -1.0;
%         ggain = .2613;
%         bgain = .0986;
        
%         %9/12/10 (10 deg fundamentals)
        rgain = -1.0;
        ggain = .3086;
        bgain = .1311;

        
    case 9 %M isoloation (Mouse)
        
        rgain = 0;
        ggain = 1;
        bgain = 0;
        
    case 10 %Luminance (S+M) isolation (Mouse)
        
        rgain = 0;
        ggain = 1;
        bgain = 1;
        
    case 11 %S  isolation (Mouse)
        
        rgain = 0;
        ggain = -1/2;
        bgain = .8934/2;
        
    case 12 %S-M isolum (Mouse)
        
        rgain = 0;
        ggain = 1;
        bgain = -.3114;
        
        
end
        
ImRGB = zeros(length(Im{1}(:,1)),length(Im{1}(1,:)),3,'uint8');  %make ImRGB uint8

if length(Im) == 1
    
    ImRdum = Im{1}*rgain;   %[-1 1]
    ImGdum = Im{1}*ggain;
    ImBdum = Im{1}*bgain;
    
else
    
    ImRdum = Im{1}*rgain;   %[-1 1]
    ImGdum = Im{2}*ggain;
    ImBdum = Im{3}*bgain;
    
end

% G = [0 0; 0 5; 5 5; 5 0;-5 0; 0 -5; -5 5; 5 -5];
% U = [.2 1; .2 6.2; 5.5 6.5; 5.2 1.2; -4.5 1; 0 -4; -5 6; 5.2 -4];
%
% trx = fitgeotrans(U,G,'affine');

% global trx
% trx2 = trx;
% pixScale = ones(size(trx2.T));
% pixScale(2,1) = 1140/912;
% pixScale(1,2) = 912/1140;
% pixScale(3,1) = xdegperpix;
% pixScale(3,2) = ydegperpix;
% 
% trx2.T = trx2.T.*pixScale
% q = imref2d(size(Im));
% [ImBdum] = imwarp(ImBdum,trx,'Outputview',q,'FillValues',0);

if ~isempty(mask) %Made 'if' to significantly reduce computation for large stimuli
    C = (128 - P.background)/128;  %[-1 1]; used so the mask looks right when background ~= 128
    ImRdum = (ImRdum + C).*mask - C; %scale w/ mask, using backround as the zero point.
    ImGdum = (ImGdum + C).*mask - C;
    ImBdum = (ImBdum + C).*mask - C;
end
ImRdum = (ImRdum+1)/2 - (.5-P.redbase);  %[0 1]
ImGdum = (ImGdum+1)/2 - (.5-P.greenbase);  
ImBdum = (ImBdum+1)/2 - (.5-P.bluebase);  

ImRGB(:,:,1) = round(ImRdum*255);  %[0 255]
ImRGB(:,:,2) = round(ImGdum*255);
ImRGB(:,:,3) = round(ImBdum*255);

