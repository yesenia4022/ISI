function synctimes = getSyncTimesIT

global ACQinfo

Nframes = ACQinfo.numberOfFrames;

synctimes = [];
for frame = 1:Nframes/3
frame
    CHsync = Load2phImage(frame,[0 0 1]);
    CHsync = CHsync{3}';
    
    stdum = getsynctimes(CHsync(:))/1000;
    
    synctimes = [synctimes; stdum(:)];

end