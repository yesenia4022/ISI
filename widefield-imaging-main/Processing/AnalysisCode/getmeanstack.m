function muStack = getmeanstack(filepath)

tf = imformats('tif');
info = feval(tf.info, filepath);
infoH = info(1).ImageDescription;
imgHeader = parseHeaderNew(infoH);
ACQinfo = imgHeader.acq;

Nslice = ACQinfo.numberOfZSlices;
dstep = ACQinfo.zStepSize;
Nframes = ACQinfo.numberOfFrames;

Nimages = length(info);

tf = imformats('tif');

framestart = 1;
framestop = Nimages;

clear CHs CHsdum
k = 0;
slce = 1;
for frame=framestart:2:framestop
    k = k+1;
    CHsdum(:,:,k) = single(feval(tf.read, filepath, frame));
    
    if k == Nframes
        k = 0;
        %stdStack{slce} = std(CHsdum,[],3);
        muStack{1}(:,:,slce) = mean(CHsdum,3);
        slce = slce+1;
    end
    
end

clear CHs CHsdum
k = 0;
slce = 1;
for frame=framestart+1:2:framestop
    k = k+1;
    CHsdum(:,:,k) = single(feval(tf.read, filepath, frame));
    
    if k == Nframes
        k = 0;
        %stdStack{slce} = std(CHsdum,[],3);
        muStack{2}(:,:,slce) = mean(CHsdum,3);
        slce = slce+1;
    end
    
end