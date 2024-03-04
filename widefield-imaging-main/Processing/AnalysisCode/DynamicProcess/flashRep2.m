function [kern popResp CoM] = flashRep2(CHs,bwCell1)

%2 uses the cell mask to create the time courses

global ACQinfo

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);
Ncell = length(celldom);

dtau = 256;
taudom = 0:dtau:15000;

rows = ACQinfo.linesPerFrame;
cols = ACQinfo.pixelsPerLine;
sp = ACQinfo.msPerLine*ACQinfo.linesPerFrame;
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

popResp = cell(1,Ncell);

for r = 1:length(CHs)
    Nframes = length(CHs{r}{1}(1,1,:));
    Tper = pepgetparam('t_period');  %stimulus period in frames
    Tper = Tper(1)*.01*100/59.55;  %60Hz frame rate

    synccourse = [];
    for i = 1:Nframes
        dum = CHs{r}{2}(:,:,i)';
        synccourse = [synccourse; dum(:)];
    end

    synctimes = getsynctimes(synccourse);

    tdom = 0:Nframes-1;
    tdom = tdom*ACQinfo.msPerLine*ACQinfo.linesPerFrame;
    fdom = linspace(0,1,length(tdom)+1);
    fdom = fdom(1:end-1)/sp*1000;

    
    for p = 1:Ncell
        
        if r == 1
            popResp{p} = cell(1,length(taudom));
        end
        
        [idcelly idcellx] = find(masklabel == celldom(p));
        idcell = find(masklabel(:) == celldom(p));
        CoM{p} = [mean(idcelly) mean(idcellx)];  %center of mass                 
        
        for z = 1:Nframes
            CHsdum = CHs{r}{1}(:,:,z);
            %tcourse(z) = (mean(CHsdum(idcell)) - mean(CHsdum(:)))/std(CHsdum(:));
            tcourse(z) = mean(CHsdum(idcell));
        end

        tcourse = zscore(tcourse);

        %tcourse = LFPfilt(tcourse',0,1000/sp,0.8,0.2)';

        fcourse = abs(fft(tcourse));
        
        tau_xy = (CoM{p}(1)-1)*ACQinfo.msPerLine + ptime*CoM{p}(2);  %approximate delay for this cell 
        tdom = (0:length(tcourse)-1)*sp;
        tdom_pix = tdom + tau_xy;   %time domain of the given cell

        % hh = zeros(1,length(tcourse));
        % hh(1:6) = ones(1,6);
        % hh = hh/sum(hh);
        % tcourse2 = ifft(fft(tcourse).*abs(fft(hh)));


        for k = 1:length(taudom)

            for j = 1:length(synctimes)
                id = find(tdom_pix>synctimes(j)+taudom(k)-dtau/2 & tdom_pix<synctimes(j)+taudom(k)+dtau/2);
                if ~isempty(id)
                    popResp{p}{k} = [popResp{p}{k} tcourse(id(1))];
                else
                    popResp{p}{k} = [popResp{p}{k} NaN];
                end

            end

        end
        

    end


end

figure
for p = 1:Ncell
    for k = 1:length(taudom)
        kern{p}(k) = nanmean(popResp{p}{k});
    end
        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
        plot(taudom,kern{p},'.-')
        %ylim([-.5 1.5])
end


