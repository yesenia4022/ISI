function Gplotretcoverage(angx,magx,angy,magy)       

%angx and angy are in degrees of visual field.  
%Important: angx goes from left to right, and angy goes from bottom to top
%to match the x_pos and y_pos domain

global bw Analyzer

screenDist = Analyzer.M.screenDist;
screenResX = Analyzer.M.xpixels/Analyzer.M.screenXcm;  %pix/cm
screenResY = Analyzer.M.ypixels/Analyzer.M.screenYcm;

%Gets position and width of bars for each condition, excluding blanks:
[xpos ypos xsize ysize] = getPosSize;

xposdom = sort(xpos);
xposdom(isnan(xposdom)) = [];
xsize_cm = (xposdom(end)-xposdom(1)+min(xsize))/screenResX;
%xsize_deg = 2*atan2(xsize_cm/2,screenDist)*180/pi;
xsize_deg = 360*xsize_cm/(2*pi*screenDist);

yposdom = sort(ypos);
yposdom(isnan(yposdom)) = [];
ysize_cm = (yposdom(end)-yposdom(1)+min(ysize))/screenResY;
%ysize_deg = 2*atan2(ysize_cm/2,screenDist)*180/pi;  %convert to deg
ysize_deg = 360*ysize_cm/(2*pi*screenDist);

%xrange = [xposdom(1)-min(xsize)/2 xposdom(end)+min(xsize)/2];
%yrange = [yposdom(1)-min(ysize)/2 yposdom(end)+min(ysize)/2];


Nx = 512; Ny = 512;

angx = (Nx-1)*angx/xsize_deg;  %0 to Nx-1
angx = round(angx)+1; %1 to Nx

angy = (Ny-1)*angy/ysize_deg;  %0 to Ny-1
angy = round(angy)+1; %1 to Ny

%Compute the moments
id = find(bw);
probX = magx(id)/sum(magx(id));
probY = magy(id)/sum(magy(id));
CoMx = sum(angx(id).*probX);
CoMy = sum(angy(id).*probY);
Varx = sum(probX.*(angx(id)-CoMx).^2);
Vary = sum(probY.*(angy(id)-CoMy).^2);

CoMx = CoMx/Nx*Analyzer.M.xpixels;
CoMy = CoMy/Ny*Analyzer.M.ypixels;
Stdx = sqrt(Varx)/Nx*xsize_deg;
Stdy = sqrt(Vary)/Ny*ysize_deg;

implot = zeros(Ny,Nx);

posID = angy(:) + (angx(:)-1)*Ny;  %vectorized location

locality = sqrt(magx + magy);

locality = locality-min(locality(:));
locality = locality/max(locality(:));
locality = locality.*bw;

id = find(~isnan(posID));
implot(posID(id)) = locality(id);

imagesc(implot), colormap gray
Xticklocs = 0:100:Nx; 
Yticklocs = 0:100:Ny;
set(gca,'Ytick',Yticklocs,'Xtick',Xticklocs);

Xticklabs = round(screenResX*xsize_cm*Xticklocs/Nx);
Yticklabs = round(screenResY*ysize_cm*Yticklocs/Ny);


set(gca,'XTickLabel',{num2str(Xticklabs(1)),num2str(Xticklabs(2)),num2str(Xticklabs(3)),num2str(Xticklabs(4)),num2str(Xticklabs(5)),num2str(Xticklabs(6))})
set(gca,'YTickLabel',{num2str(Yticklabs(1)),num2str(Yticklabs(2)),num2str(Yticklabs(3)),num2str(Yticklabs(4)),num2str(Yticklabs(5)),num2str(Yticklabs(6))})
xlabel('x position (pixels)'),ylabel('y position (pixels)')

title(['X=' num2str(round(CoMx)) ' ; Y=' num2str(round(CoMy)) ' ; Xsig=' num2str(round(Stdx)) ' deg ; Ysig=' num2str(round(Stdy)) ' deg '])