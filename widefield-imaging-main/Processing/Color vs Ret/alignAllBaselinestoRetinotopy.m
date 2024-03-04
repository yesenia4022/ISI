%function alignAllBaselinestoRetinotopy

anim = 'ry5';
dataroot = 'F:\ISIdata\';
retexp = 'u001_000';
expList = {'u001_008','u001_007','u001_006','u001_004','u001_001','u001_011'};

% anim = 'ra5';
% dataroot = 'E:\ISIdata\';
% retexp = 'u003_000'; %Kalatsky retinotopy
% %expList = {'u003_001','u001_002','u001_003','u001_004''u001_005'}; %baseline, highest to lowest
% expList = {'u002_006','u002_002','u002_003','u002_004','u002_005'};

% anim = 'rc2';
% dataroot = 'E:\ISIdata\';
% retexp = 'u002_002';
% expList = {'u002_003','u002_006','u002_007','u002_008','u002_012'}; %baseline, highest to lowest

% anim = 'rb6';
% dataroot = 'G:\ISIdata\';
% retexp = 'u003_000';
% expList = {'u003_002','u003_005','u003_006','u003_007','u003_010'}; %baseline, highest to lowest

% anim = 'ra7';
% dataroot = 'G:\ISIdata\';
% retexp = 'u001_009';
% expList = {'u001_001','u001_007','u001_003','u001_004','u001_005'}; %baseline, highest to lowest

% anim = 'rf5';
% dataroot = 'I:\ISIdata\';
% retexp = 'u000_000';
% % expList = {'u000_002','u000_005','u000_006','u000_007','u000_010'}; %baseline, highest to lowest
% expList = {'u000_001','u000_002'} % dilation experiment alone
% 

% anim = 'rf6';
% dataroot = 'I:\ISIdata\';
% retexp = 'u000_000';
% % expList = {'u000_003','u000_006','u000_007','u000_008','u000_011'}; %baseline, highest to lowest
% expList = {'u000_001','u000_003'} % dilation experiment alone
    
% anim = 'rf7';
% dataroot = 'G:\ISIdata\';
% retexp = 'u000_001';
% expList = {'u000_003','u000_004','u000_005','u000_008','u000_009'}; %baseline, highest to lowest
% 
% anim = 'rf4';
% dataroot = 'I:\ISIdata\';
% retexp = 'u000_001';
% expList = {'u000_003','u000_006','u000_008','u000_009','u000_012'}; %baseline, highest to lowest

% anim = 'rf3';
% dataroot = 'E:\ISIdata\';
% retexp = 'u000_002';
% expList = {'u000_004','u000_007','u000_008','u000_009','u000_012'}; %baseline, highest to lowest

% anim = 'rf9';
% dataroot = 'E:\ISIdata\';
% retexp = 'u000_000';
% expList = {'u000_001','u000_004','u000_005','u000_006','u000_009'}; %baseline, highest to lowest

% anim = 'rb8';
% dataroot = 'G:\ISIdata\';
% retexp = 'u002_000';
% expList = {'u002_010','u002_003','u002_005','u002_006','u002_007','u002_008'}; %baseline, highest to lowest

% anim = 'nz0';
% dataroot = 'G:\ISIdata\';
% retexp = 'u003_001';
% expList = {'u003_003','u003_014','u003_015','u003_016','u003_022'}; %baseline, highest to lowest

% anim = 'nz0';
% dataroot = 'G:\ISIdata\';
% retexp = 'u003_001';
% expList = {'u003_002','u003_014','u003_015','u003_016','u003_022'}; %baseline, highest to lowest

% anim = 'rc3';
% dataroot = 'E:\ISIdata\';
% retexp = 'u006_000';
% % expList = {'u006_011','u006_010','u006_004','u006_005','u006_006','u006_009'}; %baseline, highest v2 to lowest
% % expList = {'u006_014','u006_015','u006_017','u006_018','u006_019','u006_020'}; %2ND complete run v3
% expList = {'u006_014','u006_010','u006_017','u006_018','u006_006','u006_020'}; %Combo of v2 and v3

% anim = 'rc4';
% dataroot = 'g:\ISIdata\';
% retexp = 'u005_001';
% expList = {'u005_014','u005_003','u005_004','u005_005','u005_006','u005_009'};
% 
% anim = 'rk3';
% dataroot = 'e:\ISIdata\';
% retexp = 'u003_020';
% expList = {'u003_002','u003_003','u003_006','u003_007','u003_010','u003_011'};

% anim = 'rk4';
% dataroot = 'f:\ISIdata\';
% retexp = 'u002_000';
% expList = {'u002_016','u002_002','u002_005','u002_014','u002_008','u002_011'};

% anim = 'rk7';
% dataroot = 'E:\ISIdata\';
% retexp = 'u002_000';
% expList = {'u002_001','u002_002','u002_003','u002_004','u002_005','u002_006'};

% anim = 'rl5';
% dataroot = 'h:\ISIdata\';
% retexp = 'u004_010';
% expList = {'u005_001','u004_000','u004_002','u004_003','u004_004','u004_005'};

% anim = 'rl5'; %FOR S/T Experiment
% dataroot = 'h:\ISIdata\';
% retexp = 'u004_010';
% expList = {'u005_003'};


% anim = 'rl7'; % first set of baselines
% dataroot = 'E:\ISIdata\';
% retexp = 'u003_001';
% expList = {'u003_002','u003_003','u003_004','u003_005','u003_006','u003_007'};% experiments in order, but rerun of B5 is better. Might need to use exp 8 for B5

% anim = 'rl7'; % second set of baselines
% dataroot = 'E:\ISIdata\';
% retexp = 'u003_001';
% expList = {'u003_009','u003_011','u003_012','u003_013','u003_014','u003_015'};

% 
% anim = 'rl8'; % 
% dataroot = 'H:\ISIdata\';
% retexp = 'u002_002';
% expList = {'u002_004','u002_005','u002_006','u002_008','u002_009','u002_010'};

% anim = 'rm8'; %
% dataroot = 'i:\ISIdata\';
% retexp = 'u003_000';
% expList = {'u003_004','u002_005','u002_006','u002_008','u002_009','u002_010'};

% anim = 'rm2'; % linear set of baselines
% dataroot = 'i:\ISIdata\';
% retexp = 'u002_000';
% expList = {'u002_008','u002_001','u002_002','u002_003','u002_006','u002_007'};

% anim = 'rm2'; % log set of baselines v1
% dataroot = 'i:\ISIdata\';
% retexp = 'u002_000';
% expList = {'u002_008','u002_001','u002_011','u002_012','u002_013','u002_014'};
%
% anim = 'rn1'; % log set of baselines v1
% dataroot = 'F:\ISIdata\';
% retexp = 'u003_000';
% expList = {'u003_007','u003_001','u003_002','u003_003','u003_004','u003_006'};

% anim = 'rm8'; % log set of baselines v1
% dataroot = 'F:\ISIdata\';
% retexp = 'u003_000';
% expList = {'u002_007','u002_002','u002_003','u002_004','u002_005','u002_006'};

% anim = 'rm9'; % log set of baselines missing B6, problems with B1
% dataroot = 'H:\ISIdata\';
% retexp = 'u001_000';
% expList = {'u002_003','u002_005','u002_007','u002_008','u002_012'};

% anim = 'rj9'; % log set of baselines v1
% dataroot = 'h:\ISIdata\';
% retexp = 'u00_000';
% expList = {'u002_007','u002_002','u002_003','u002_004','u002_005','u002_006'};

% anim = 'ra6'; % log set of baselines v1
% dataroot = 'h:\ISIdata\';
% retexp = 'u00_000';
% expList = {'u002_007','u002_002','u002_003','u002_004','u002_005','u002_006'};

% anim = 'rm0'; % log set of baselines using U0 retinotopy  
% dataroot = 'h:\ISIdata\';
% retexp = 'u000_000';
% expList = {'u002_007','u002_002','u002_003','u002_004','u002_005','u002_006'};

% anim = 'rm0'; % linear set of baselines using U0 retinotopy 
% dataroot = 'h:\ISIdata\';
% retexp = 'u000_000';
% expList = {'u002_007','u002_002','u002_008','u002_009','u002_010','u002_011'};

% anim = 'rm0'; % log set of baselines using U2 retinotopy  
% dataroot = 'h:\ISIdata\'; %gtrx_rm0 v4
% retexp = 'u002_013';
% expList = {'u002_007','u002_012','u002_003','u002_004','u002_005','u002_006'};

% anim = 'rm0'; % linear set of baselines using U2 retinotopy 
% dataroot = 'h:\ISIdata\'; %gtrx_rm0 v5
% retexp = 'u002_013';
% expList = {'u002_007','u002_012','u002_008','u002_009','u002_010','u002_011'};

% anim = 'sb6'; % linear set of baselines using U2 retinotopy 
% dataroot = 'E:\ISIdata\'; %gtrx_rm0 v5
% retexp = 'u002_000';
% expList = {'u002_014','u002_003','u002_006','u002_007','u002_008','u002_009'};

% anim = 'ry5'; % linear set of baselines using U2 retinotopy 
% dataroot = 'E:\ISIdata\'; %gtrx_rm0 v5
% retexp = 'u001_000';
% expList = {'u001_011','u001_001','u001_004','u001_006','u001_007','u001_008'};

% anim = 'sf6'; % linear set of baselines using U2 retinotopy 
% dataroot = 'f:\ISIdata\'; %gtrx_rm0 v5
% retexp = 'u001_000';
% expList = {'u001_002','u001_005','u001_006','u001_007','u001_008'};

% anim = 'se3'; % linear set of baselines using U2 retinotopy 
% dataroot = 'E:\ISIdata\'; %gtrx_rm0 v5
% retexp = 'u001_002';
% expList = {'u001_010','u001_003','u001_006','u001_007','u001_008','u001_009'};

% anim = 'se4'; % linear set of baselines using U2 retinotopy 
% dataroot = 'E:\ISIdata\'; %gtrx_rm0 v5
% retexp = 'u001_001';
% expList = {'u001_013','u001_002','u001_009','u001_010','u001_011', 'u001_012'};


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imloc = [dataroot anim '\' retexp '\' retexp '_000_f50.mat'];
S = load(imloc);  %Returns the contents in the .mat under the structure S
im_output = double(S.im);

im_output = im_output-prctile(im_output(:),1);
im_output(find(im_output<0)) = 0;
im_output = im_output/prctile(im_output(:),99.9);
im_output(find(im_output>1)) = 1;

Nbase = length(expList);

for bL = 1:Nbase
        
    disp(['Aligning baseline experiment ' num2str(bL) ' to retinotopy'])
    
    imloc = [dataroot anim '\' expList{bL} '\' expList{bL} '_000_f50.mat'];
    S = load(imloc);  %Returns the contents in the .mat under the structure S
    im_input = double(S.im);
    
    im_input = im_input-prctile(im_input(:),1);
    im_input(find(im_input<0)) = 0;
    im_input = im_input/prctile(im_input(:),99.9);
    im_input(find(im_input>1)) = 1;
     
    [movingPoints,fixedPoints] = cpselect(im_input,im_output,'Wait',true);
    
    trx{bL} = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
       
    im_in{bL} = im_input;
end

%%
path = 'c:\baselineTrx';
uisave('trx',path)


%%

[xp yp] = meshgrid(1:20:size(im_output,2),1:20:size(im_output,1));
figure
subplot(1,Nbase+1,1)
imagesc(im_output), colormap gray,axis image
hold on, plot(xp(:),yp(:),'.r')
title('retinotopy experiment')
for bL = 1:Nbase
    Ref_Im = imref2d(size(im_output));
    dum = imwarp(im_in{bL},trx{2},'OutputView',Ref_Im);
    subplot(1,Nbase+1,bL+1)
    imagesc(dum), colormap gray,  axis image
    
    hold on, plot(xp(:),yp(:),'.r')
end




% 
% 
% Ref_Im = imref2d(size(im_input));
% 
% imstate.fmaps{1} = imwarp(imstate.fmaps{1},trx,'OutputView',Ref_Im);
% imstate.fmaps{2} = imwarp(imstate.fmaps{2},trx,'OutputView',Ref_Im);
% imstate.sigMag = imwarp(imstate.sigMag,trx,'OutputView',Ref_Im);
% imstate.bw = imwarp(imstate.bw,trx,'OutputView',Ref_Im);
% imstate.imfunc = imwarp(imstate.imfunc,trx,'OutputView',Ref_Im);
% imstate.imanatFunc = imwarp(imstate.imanatFunc,trx,'OutputView',Ref_Im);
% imstate.mag = imwarp(imstate.mag,trx,'OutputView',Ref_Im);
% imstate.imanat = imwarp(imstate.imanat,trx,'OutputView',Ref_Im);
% 
% if ~isempty(imstate.areaBounds)
%     imstate.areaBounds = floor(imwarp(imstate.areaBounds,trx,'OutputView',Ref_Im));
%     SE = strel('disk',1);
%     imstate.areaBounds = imopen(imstate.areaBounds,SE); %Trx tends to fuse areas
% end
% 
% Im_registered = imwarp(im_input,trx,'OutputView',Ref_Im);
% 
% figure,
% subplot(1,3,1), imagesc(im_input), title('unregistered input image'), axis image
% for i =1:length(fixedPoints(:,1))
%     hold on,
%     plot(fixedPoints(i,1),fixedPoints(i,2),'xr')
% end
% 
% 
% subplot(1,3,2), imagesc(Im_registered), title('registered input image'), axis image
% for i =1:length(fixedPoints(:,1))
%     hold on,
%     plot(fixedPoints(i,1),fixedPoints(i,2),'xr')
% end
% 
% subplot(1,3,3), imagesc(im_output), title('output image'), axis image
% for i =1:length(fixedPoints(:,1))
%     hold on,
%     plot(fixedPoints(i,1),fixedPoints(i,2),'xr')
% end
% colormap gray
