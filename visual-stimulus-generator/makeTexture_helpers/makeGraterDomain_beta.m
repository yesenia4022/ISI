function [sdom, tdom, x_ecc, y_ecc] = makeGraterDomain_beta(xN,yN,ori,s_freq,t_period,altazimuth,imtrx)

%Ian Nauhaus

global Mstate screenNum

P = getParamStruct;

affinebit = 1;

screenRes = Screen('Resolution',screenNum);

pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

%% First identify the rectangular patch on the screen, along with the pixel
%%coordinate system within the patch

xpos = P.x_pos;
ypos = P.y_pos;

%Get upper and lower limits of the patch in screen pixel coordinates
xlim = [xpos - xN/2*P.x_zoom    xpos + xN/2*P.x_zoom]; 
ylim = [ypos - yN/2*P.y_zoom    ypos + yN/2*P.y_zoom];

x_ecc = single(linspace(xlim(1),xlim(2),xN));  %pixel domain on the screen 
y_ecc = single(linspace(ylim(1),ylim(2),yN));  

[x_ecc y_ecc] = meshgrid(x_ecc,y_ecc);  %pix

if affinebit
    x_eccdum = x_ecc; y_eccdum = y_ecc; oneMat = ones(size(y_ecc));
    x_ecc = imtrx(1,1)*x_eccdum + imtrx(2,1)*y_eccdum + imtrx(3,1)*y_eccdum.*x_eccdum + imtrx(4,1)*oneMat;
    y_ecc = imtrx(1,2)*x_eccdum + imtrx(2,2)*y_eccdum + imtrx(3,2)*y_eccdum.*x_eccdum + imtrx(4,2)*oneMat;
else
    x_eccdum = x_ecc; y_eccdum = y_ecc; oneMat = ones(size(y_ecc));
    x_ecc = imtrx(1,1)*x_eccdum + imtrx(2,1)*y_eccdum + imtrx(3,1)*oneMat;
    y_ecc = imtrx(1,2)*x_eccdum + imtrx(2,2)*y_eccdum + imtrx(3,2)*oneMat;
end

x_ecc = x_ecc-xpos;
y_ecc = y_ecc-ypos;

%Scale to convert pixel domain to cm or deg domain
switch altazimuth

    case 'none'  %For this case we assume that the screen is curved
        
        xdegperpix = P.x_size/(xN*P.x_zoom);
        ydegperpix = P.y_size/(yN*P.y_zoom);
        
        x_ecc = x_ecc*xdegperpix; %This scaling is the same for affine or not
        y_ecc = y_ecc*ydegperpix; %deg
       
    case 'altitude'
        
        %%%Get the xy domain
        
        x_ecc = x_ecc/pixpercmX; %cm
        y_ecc = y_ecc/pixpercmY;
        
        
    case  'azimuth'
        
        x_ecc = x_ecc/pixpercmX; %cm
        y_ecc = y_ecc/pixpercmY;
             
end

%Change location of perpendicular bisector relative to stimulus center
x_ecc = x_ecc-P.dx_perpbis;
y_ecc = y_ecc-P.dy_perpbis;
%%%%%%%%%%%%%%%%%%%%%%%


switch altazimuth
    
    case 'none'
        
        sdom = x_ecc*cos(ori*pi/180) - y_ecc*sin(ori*pi/180);    %deg
        
    case 'altitude'
        
        %Apply "tilt" to y/z dimensions: rotation around x axis
        z_ecc = Mstate.screenDist*ones(size(x_ecc));  %dimension perpendicular to screen
        y_eccT = y_ecc*cos(P.tilt_alt*pi/180) - z_ecc*sin(P.tilt_alt*pi/180);
        z_eccT = y_ecc*sin(P.tilt_alt*pi/180) + z_ecc*cos(P.tilt_alt*pi/180);       
             
        %Apply "tilt direction", i.e. rotation around y axis   
        x_eccR = x_ecc*cos(P.tilt_az*pi/180) - z_eccT*sin(P.tilt_az*pi/180);
        z_eccR = x_ecc*sin(P.tilt_az*pi/180) + z_eccT*cos(P.tilt_az*pi/180); 
        
        %Apply "orientation" to the x/y dimensions: rotation around z axis
        x_eccO = x_eccR*cos(ori*pi/180) - y_eccT*sin(ori*pi/180); 
        y_eccO = x_eccR*sin(ori*pi/180) + y_eccT*cos(ori*pi/180); 
        
        xdum = y_eccO; %If I don't swap them then the orientation is wrong by 90deg.
        ydum = x_eccO;
                
        sdom = asin(ydum./sqrt(xdum.^2 + ydum.^2 + z_eccR.^2))*180/pi;

    case 'azimuth' %The projection of azimuth onto a plane is the same as a cylinder on a plane       
                
        %Apply "tilt" to y/z dimensions: rotation around x axis
        z_ecc = Mstate.screenDist*ones(size(x_ecc));  %dimension perpendicular to screen
        y_eccT = y_ecc*cos(P.tilt_alt*pi/180) - z_ecc*sin(P.tilt_alt*pi/180);
        z_eccT = y_ecc*sin(P.tilt_alt*pi/180) + z_ecc*cos(P.tilt_alt*pi/180);       
             
        %Apply "tilt direction", i.e. rotation around y axis   
        x_eccR = x_ecc*cos(P.tilt_az*pi/180) - z_eccT*sin(P.tilt_az*pi/180);
        z_eccR = x_ecc*sin(P.tilt_az*pi/180) + z_eccT*cos(P.tilt_az*pi/180); 
        
        %Apply "orientation" to the x/y dimensions: rotation around z axis
        x_eccO = x_eccR*cos(ori*pi/180) - y_eccT*sin(ori*pi/180); 
        
        sdom = atan(x_eccO./z_eccR)*180/pi; %deg
        
end

sdom = sdom*s_freq*2*pi; %radians

if t_period == inf   
    tdom = 0;
else
    tdom = single(linspace(0,2*pi,t_period+1));
    tdom = tdom(1:end-1);
end
    
    
    

    