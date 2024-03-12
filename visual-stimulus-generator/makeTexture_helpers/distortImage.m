function ImD = distortImage(Im,P)

global Mstate

%P.x_size = 160; P.y_size = 160;


xcm = 2*Mstate.screenDist*tan(P.x_size/2*pi/180);  %width in cm
ycm = 2*Mstate.screenDist*tan(P.y_size/2*pi/180);  %height in cm

xdom = linspace(0,xcm,size(Im,2));
ydom = linspace(0,ycm,size(Im,1));

[x y] = meshgrid(xdom,ydom);

% x = x-mean(xdom(:))-P.perpBisX;
% y = y-mean(ydom(:))-P.perpBisY;
x = x-mean(xdom(:));
y = y-mean(ydom(:));


%Create map of the spherical coordinates within presentation window: [theta(x,y) phi(x,y)]
ang = atan2(y,x);
rho = atan2(sqrt(x.^2 + y.^2) , Mstate.screenDist);
xp = rho.*cos(ang)*180/pi;
yp = rho.*sin(ang)*180/pi;

%Create map of the same spherical domain, but on the image (n.b. must be a grid to run interp2)

xdom = linspace(0,P.x_size,size(Im,2));
ydom = linspace(0,P.y_size,size(Im,1));
[x y] = meshgrid(xdom-mean(xdom),ydom-mean(ydom));
ImD = interp2(x,y,Im,xp,yp);


%figure,imagesc(ImD)




