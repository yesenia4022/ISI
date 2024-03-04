function [M_Contrast S_Contrast] = getConeContrast(gdom,udom)


%% Get Spectral measurements and sensitivity functions

root = 'C:\Stimulator_master\Calibration\UVprojector\5_8_16\';

f = [root 'spectrum_UV128'];  
load(f,'I','dom')
UV = I;

%In 5_8_16 calibration, I set the current of the green during calibration to 65
root = 'C:\Stimulator_master\Calibration\RGBprojector\5_8_16\';

f = [root 'spectrum_green128'];
load(f,'I','dom')
%%%%%%%%%%%N.B>>>>
green = I;  

% figure,plot(dom,green)
% hold on,plot(dom,UV)

% Get sensitivity functions

[rod,Mcone,Scone,lambda] = photoreceptor_templates(500,508,360);
id = find(lambda>=dom(1) & lambda<=dom(end));
lambda = lambda(id);
Scone = Scone(id);
Mcone = Mcone(id);

UV   = interp1(dom,UV,lambda,'spline');
green = interp1(dom,green,lambda,'spline');


%%

R = 1;

Ti = [21.1055 -38.8909; -0.0297 25.6989]; %'Middle' Calibrated 5/8/16
%Ti = [35.6516 -65.9533; -0.1890 38.9525]; %'Right' Calibrated 5/8/16

contrastNorm = 1;
if contrastNorm
    T=inv(Ti); %Bring it back to the noninverted matrix
    alpha = sum(T(1,:))/sum(T(2,:)); %Ratio of the cone "responses" at the gray point.
    Ti(:,2) = Ti(:,2)/alpha; %Correct for the different gray points of the two cones
end

%What I was doing
% Ti(1,:) = Ti(1,:)/norm(Ti(1,:));
% Ti(2,:) = Ti(2,:)/norm(Ti(2,:));

%What I should be doing
normer = max([norm(Ti(1,:)) norm(Ti(2,:))]);
Ti(1,:) = Ti(1,:)/normer;
Ti(2,:) = Ti(2,:)/normer;


i = 0;

clear TotalContrast G_UV_Contrast ang

for u = 1:length(udom)
    for g = 1:length(gdom)
        
        %     M = R*cos(thetadom(i)*pi/180);
        %     S = R*sin(thetadom(i)*pi/180);
        %     Liso = Ti*[M S]';
        
        Liso = [gdom(g) udom(u)];
        
        
        %This next stuff replicates (line-by-line) stimulus generation code (ImtoRGB.m)
        
        Glin = linspace(0,2,256) - 1;
        Ulin = linspace(0,2,256) - 1; %[-1 1]
        
        Glin = Glin*Liso(1);
        Ulin = Ulin*Liso(2); %[-1 1]
        
        Glin = (Glin + 1)/2;
        Ulin = (Ulin + 1)/2;  %[0 1]
        
        %N.B. Rounding will create slightly imperfect isolation
        % Ulin = round(Ulin*255);  %[0 255]
        % Glin = round(Glin*255);
        
        %Estimate min/mid/max spectra using scaling coefficients
        %Spectra were computed at the maximum buffer value
        baselinei = 0;
        SPECmin = [green(:) UV(:)]*[Glin(1) Ulin(1)]' - baselinei; %spectrum at trough (assumes baseline was subtracted up top)
        SPECmax = [green(:) UV(:)]*[Glin(end) Ulin(end)]' -baselinei; %spectrum at peak
        
        
        Rmax = [Mcone(:) Scone(:)]'*SPECmax;
        Rmin = [Mcone(:) Scone(:)]'*SPECmin;
        Rmean = (Rmax+Rmin)/2;
        dum = (Rmax-Rmin)./(2*Rmean);
        M_Contrast(u,g) = dum(1);
        S_Contrast(u,g) = dum(2);
        ang(u,g) = atan2(S_Contrast(u,g),M_Contrast(u,g))*180/pi;
        
    end
end


