function CHs = getStack(filepath,chvec)

%e.g. CHs = getStack('C:\2p_data\fd9\fd9fast002.tif',[1 1 0 1]);%

tf = imformats('tif');
info = feval(tf.info, filepath);
infoH = info(1).ImageDescription;
imgHeader = parseHeaderNew(infoH);
ACQinfo = imgHeader.acq;

Nframes = ACQinfo.numberOfFrames;

acqId = [ACQinfo.acquiringChannel1 ACQinfo.acquiringChannel2 ACQinfo.acquiringChannel3 ACQinfo.acquiringChannel4]; 

nZ = ACQinfo.numberOfZSlices;

dFrame = ACQinfo.numberOfChannelsAcquire;

dim = [ACQinfo.linesPerFrame ACQinfo.pixelsPerLine];


if chvec(1)   
    for z = 1:nZ  %loop through each slice
        fstart = (z-1)*dFrame*Nframes+1;        
        fstop = z*dFrame*Nframes;
        dim2 = [dim length(fstart:dFrame:fstop)];
        CHs{1}{z} = zeros(dim2,'single');
        k = 1;
        for frame = fstart:dFrame:fstop
            A = feval(tf.read, filepath, frame);
            CHs{1}{z}(:,:,k) = single(A);
            k = k+1;
        end
    end
end

if chvec(2)
    for z = 1:nZ  %loop through each slice
        fstart = (z-1)*dFrame*Nframes+1 + acqId(1);        
        fstop = z*dFrame*Nframes;
        dim2 = [dim length(fstart:dFrame:fstop)];
        CHs{2}{z} = zeros(dim2,'single');
        k = 1;
        for frame = fstart:dFrame:fstop
            A = feval(tf.read, filepath, frame);
            CHs{2}{z}(:,:,k) = single(A);
            k = k+1;
        end
    end
end

if chvec(3)
    for z = 1:nZ  %loop through each slice
        fstart = (z-1)*dFrame*Nframes+1 + sum(acqId(1:2));        
        fstop = z*dFrame*Nframes;
        dim2 = [dim length(fstart:dFrame:fstop)];
        CHs{3}{z} = zeros(dim2,'single');
        k = 1;
        for frame = fstart:dFrame:fstop
            A = feval(tf.read, filepath, frame);
            CHs{3}{z}(:,:,k) = single(A);
            k = k+1;
        end
    end
end

if chvec(4)
    for z = 1:nZ  %loop through each slice
        fstart = (z-1)*dFrame*Nframes+1 + sum(acqId(1:3));        
        fstop = z*dFrame*Nframes;
        dim2 = [dim length(fstart:dFrame:fstop)];
        CHs{4}{z} = zeros(dim2,'single');
        k = 1;
        for frame = fstart:dFrame:fstop
            A = feval(tf.read, filepath, frame);
            CHs{4}{z}(:,:,k) = single(A);
            k = k+1;
        end
    end
end
