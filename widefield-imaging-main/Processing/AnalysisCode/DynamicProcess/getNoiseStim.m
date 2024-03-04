function Im = getNoiseStim(cond)


global Analyzer

pixpercmX = Analyzer.M.xpixels/Analyzer.M.screenXcm;
pixpercmY = Analyzer.M.ypixels/Analyzer.M.screenYcm;

xcm = 2*pi*Analyzer.M.screenDist*getparam('x_size')/360;  %stimulus width in cm
xN = round(xcm*pixpercmX);  %stimulus width in pixels
ycm = 2*pi*Analyzer.M.screenDist*getparam('y_size')/360;   %stimulus height in cm
yN = round(ycm*pixpercmY);  %stimulus height in pixels

if strcmp(Analyzer.M.monitor,'CRT')
    refrate = 100;
elseif strcmp(Analyzer.M.monitor,'LCD')
    refrate = 60;
end

N_Im = round(getparam('stim_time')*refrate/getparam('h_per')); %number of images to present
xN = round(xN/getparam('x_zoom'));  %Downsample for the zoom
yN = round(yN/getparam('y_zoom'));

for i = 1:length(Analyzer.loops.conds{cond}.val);
    if strcmp(Analyzer.loops.conds{cond}.symbol{i},'rseed')
        rseed = Analyzer.loops.conds{cond}.val{i};
    end
end


RandStream.setDefaultStream(RandStream('mt19937ar','seed',rseed));
Im = rand(yN,xN,N_Im);  
Im = round(Im);