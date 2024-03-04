
global imfile expt

anim = 'XX0';
expt = 'u000_002';
dfile = ['c:\2photon\neurodata\' anim '\'];

pepsetdatadirectory(dfile)
pepload(expt)

imfile = [dfile expt]; 
setacqinfo(imfile);

[Data1 Data2] = F0images([0 3000]);

Omap1 = OrimapGen(Data1);
Omap2 = OrimapGen(Data2);

figure,imagesc(angle(Omap1)),colormap hsv
figure,imagesc(angle(Omap2)),colormap hsv