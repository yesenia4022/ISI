function [kern] = flashRepMovie(CHs)

%2 uses the cell mask to create the time courses

global ACQinfo

dtau = 256;
taudom = -6000:dtau:6000;

rows = ACQinfo.linesPerFrame;
cols = ACQinfo.pixelsPerLine;
imElem = rows*cols;

sp = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)       

dumIm = zeros(rows,cols);
kern = zeros(rows,cols,length(taudom));
counter = zeros(rows,cols,length(taudom));

for r = 1:length(CHs)
    Nframes = length(CHs{r}{1}(1,1,:));

    VCHs = zeros(numel(CHs{r}{1}),1);
    for z = 1:Nframes
        CHsdum = CHs{r}{1}(:,:,z)';
        %tcourse(z) = (mean(CHsdum(idcell)) - mean(CHsdum(:)))/std(CHsdum(:));
        ran = [(z-1)*imElem+1 z*imElem];
        VCHs(ran(1):ran(2)) = CHsdum(:);
    end

    Tper = pepgetparam('t_period');  %stimulus period in frames
    Tper = Tper(1)*.01*100/59.55;  %60Hz frame rate

    synccourse = [];
    for i = 1:Nframes
        dum = CHs{r}{2}(:,:,i)';
        synccourse = [synccourse; dum(:)];
    end

    synctimes = getsynctimes(synccourse);

    VCHs = zscore(VCHs);
    
    tdom = (0:length(VCHs)-1)*ptime;

    for j = 1:length(synctimes)
        for k = 1:length(taudom)

            id = find(tdom>synctimes(j)+taudom(k)-dtau/2 & tdom<synctimes(j)+taudom(k)+dtau/2);

            idf = id - imElem*(ceil(id/(imElem)) - 1);  %vector index within the frame

            dumIm = dumIm*0;
            dumIm(idf) = VCHs(id);
            kern(:,:,k) = kern(:,:,k) + dumIm';
            dumIm(idf) = 1;
            counter(:,:,k) = counter(:,:,k) + dumIm';

        end

    end


end

kern = kern./counter;



