function CH = LoadTiff2Mat(filepath)



tf = imformats('tif');

info = feval(tf.info, filepath);

Nimages = length(info);

Nch = ACQinfo.numberOfChannelsAcquire;

if frame>Nimages
    frame = Nimages;
end

if chvec(1) == 1
    CHs{1} = single(feval(tf.read, filepath, frame*Nch-(Nch-1)));
end