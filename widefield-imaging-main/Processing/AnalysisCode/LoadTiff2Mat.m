function [CHs ACQinfo] = LoadTiff2Mat(filepath,chvec)

%filepath = 'Z:\MouseDuras\EyalAD001.tif';
%filepath = 'C:\2p_data\zc4\stack000.tif';

tf = imformats('tif');
info = feval(tf.info, filepath);
infoH = info(1).ImageDescription;
imgHeader = parseHeaderNew(infoH);
ACQinfo = imgHeader.acq;

tf = imformats('tif');
info = feval(tf.info, filepath);

acqvec = [ACQinfo.acquiringChannel1 ACQinfo.acquiringChannel2 ACQinfo.acquiringChannel3 ACQinfo.acquiringChannel4];

Nimages = length(info);

fstop = Nimages;

dFrame = ACQinfo.numberOfChannelsAcquire;

dim = [ACQinfo.linesPerFrame ACQinfo.pixelsPerLine ACQinfo.numberOfFrames];

cloop = 1;
for ch = 1:4  %Loop through all 4 analog channels

    if chvec(ch) == 1
        fstart = sum(acqvec(1:ch));
        CHs{cloop} = zeros(dim,'single');
        k = 1;
        for frame=fstart:dFrame:fstop
            A = feval(tf.read, filepath, frame);
            CHs{cloop}(:,:,k) = single(A);
            k = k+1;
        end
        cloop = cloop + 1;
    end

end