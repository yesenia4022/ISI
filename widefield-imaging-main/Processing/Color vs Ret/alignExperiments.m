function alignExperiments(anim,tempexp,trxexp)

spath = 'c:\ISIgeoTrx\'; %Root location of where to save transformation structure
dataroot = 'E:\ISIdata\'; %Root location of data

%tempexp is the experiment that does not get transformed
%trxexp is the experiment that gets transformed

% anim = 'rk7';
% tempexp = 'u002_000'; %Experiment that gets aligned "to"
% % trxexp = 'u002_007'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u002_009'; %B2 sfreq
% % trxexp = 'u002_008'; %B5 tperiod
% trxexp = 'u002_010'; %B2 tperiod

% anim = 'rm1';
% tempexp = 'u001_000'; %Experiment that gets aligned "to"
% %trxexp = 'u001_002'; %B5 sfreq; Experiment that gets transformed
% trxexp = 'u001_004'; %B2 sfreq

% anim = 'rm8';
% tempexp = 'u003_000'; %Experiment that gets aligned "to"
% % trxexp = 'u003_005'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u003_009'; %B1 sfreq
% % trxexp = 'u003_004'; %B5 tperiod
% trxexp = 'u003_010'; %B1 tperiod

% anim = 'rm0';
% tempexp = 'u001_001'; %Experiment that gets aligned "to"
% % trxexp = 'u001_003'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u001_006'; %B1 sfreq
% % trxexp = 'u001_004'; %B5 tperiod
% trxexp = 'u001_008'; %B1 tperiod

% anim = 'rk4';
% tempexp = 'u003_007'; %Experiment that gets aligned "to"
% % trxexp = 'u003_002'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u003_005'; %B2 sfreq
% % trxexp = 'u003_003'; %B5 tperiod
% % trxexp = 'u003_006'; %B2 tperiod

% anim = 'rk6';
% tempexp = 'u001_002'; %Experiment that gets aligned "to"
% %  trxexp = 'u001_009'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u001_013'; %B2 sfreq
% % trxexp = 'u001_010'; %B5 tperiod
% % trxexp = 'u001_014'; %B2 tperiod

% anim = 'rl0';
% tempexp = 'u003_002'; %Experiment that gets aligned "to"
% %   trxexp = 'u003_004'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u003_011'; %B2 sfreq
% % trxexp = 'u003_005'; %B5 tperiod
% trxexp = 'u003_012'; %B2 tperiod

% anim = 'rk5';
% tempexp = 'u003_000'; %Experiment that gets aligned "to"
% %   trxexp = 'u003_002'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u003_009'; %B2 sfreq
% % trxexp = 'u003_003'; %B5 tperiod
% trxexp = 'u003_010'; %B2 tperiod

% anim = 'rl7';
% tempexp = 'u005_000'; %Experiment that gets aligned "to"
% %   trxexp = 'u005_002'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u005_005'; %B1 sfreq [0 90]
% % trxexp = 'u005_006'; %B1 sfreq [45 135]

% anim = 'rl7';
% tempexp = 'u006_000'; %Experiment that gets aligned "to"
% % trxexp = 'u006_002'; %B5 tperiod [0 90] transformed
% % trxexp = 'u003_003'; %B5 tperiod [45 135]
% % trxexp = 'u003_005'; %B1 tperiod [0 90]
% trxexp = 'u003_006'; %B1 tperiod [45 135]

% anim = 'rl5';
% tempexp = 'u005_000'; %Experiment that gets aligned "to"
% %   trxexp = 'u005_003'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u005_011'; %B1 sfreq
% % trxexp = 'u005_004'; %B5 tperiod
% trxexp = 'u005_012'; %B tperiod

% anim = 'rl5';
% tempexp = 'u006_001'; %Experiment that gets aligned "to"
% %   trxexp = 'u006_003'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u006_006'; %B1 sfreq
% % trxexp = 'u006_004'; %B5 tperiod
% % trxexp = 'u006_007'; %B1 tperiod

% anim = 'rn4';
% tempexp = 'u004_001'; %Experiment that gets aligned "to"
% %   trxexp = 'u004_004'; %B5 sfreq; Experiment that gets transformed
% trxexp = 'u004_007'; %B1 sfreq
% 
% anim = 'rk4';
% tempexp = 'u002_000'; %Experiment that gets aligned "to"
% % trxexp = 'u002_003'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u002_009'; %B2 sfreq
% % trxexp = 'u002_004'; %B5 tperiod
% % trxexp = 'u002_010'; %B2 tperiod

% anim = 'rl1';
% tempexp = 'u003_000'; %Experiment that gets aligned "to"
% % trxexp = 'u003_003'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u003_009'; %B2 sfreq
% % trxexp = 'u003_005'; %B5 tperiod

% anim = 'rm4';
% tempexp = 'u001_001'; %Experiment that gets aligned "to"
% % trxexp = 'u001_003'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u001_008'; %B2 sfreq
% trxexp = 'u001_007'; %B2 tperiod

% anim = 'nz0';
% tempexp = 'u003_001'; %Experiment that gets aligned "to"
% % trxexp = 'u003_012'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u003_020'; %B2 sfreq
% % trxexp = 'u003_013'; %B5 tperiod
% % trxexp = 'u003_021'; %B2 tperiod

% anim = 'rc2';
% tempexp = 'u002_002'; %Experiment that gets aligned "to"
% % trxexp = 'u002_004'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u002_009'; %B2 sfreq
% % trxexp = 'u002_005'; %B5 tperiod
% trxexp = 'u002_010'; %B2 tperiod

%anim = 'rb6';
%tempexp = 'u003_000'; %Experiment that gets aligned "to"
% % trxexp = 'u003_003'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u003_008'; %B2 sfreq
% % trxexp = 'u003_004'; %B5 tperiod
% % trxexp = 'u003_009'; %B2 tperiod

% anim = 'rf5';
% tempexp = 'u000_000'; %Experiment that gets aligned "to"
% % trxexp = 'u000_003'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u000_008'; %B2 sfreq
% % trxexp = 'u000_004'; %B5 tperiod
% % trxexp = 'u000_009'; %B2 tperiod

% anim = 'rf6';
% tempexp = 'u000_000'; %Experiment that gets aligned "to"
% % trxexp = 'u000_004'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u000_009'; %B2 sfreq
% % trxexp = 'u000_005'; %B5 tperiod
% % trxexp = 'u000_010'; %B2 tperiod

% anim = 'rf4';
% tempexp = 'u000_001'; %Experiment that gets aligned "to"
% % trxexp = 'u000_004'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u000_010'; %B2 sfreq
% % trxexp = 'u000_005'; %B5 tperiod
% % trxexp = 'u000_011'; %B2 tperiod

% anim = 'rf3';
% tempexp = 'u000_002'; %Experiment that gets aligned "to"
% % trxexp = 'u000_005'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u000_010'; %B2 sfreq
% % trxexp = 'u000_006'; %B5 tperiod
% % trxexp = 'u000_011'; %B2 tperiod

% anim = 'rf9';
% tempexp = 'u000_000'; %Experiment that gets aligned "to"
% % trxexp = 'u000_002'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u000_007'; %B2 sfreq
% % trxexp = 'u000_003'; %B5 tperiod
% % trxexp = 'u000_008'; %B2 tperiod

% anim = 'rb8';
% tempexp = 'u002_000'; %Experiment that gets aligned "to"
% % trxexp = 'u002_012'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u002_014'; %B2 sfreq
% % trxexp = 'u002_013'; %B5 tperiod
% % trxexp = 'u002_015'; %B2 tperiod

% anim = 'rc3';
% tempexp = 'u006_000'; %Experiment that gets aligned "to"
% % trxexp = 'u006_002'; %B5 sfreq; Experiment that gets transformed
% % trxexp = 'u006_007'; %B2 sfreq
% % trxexp = 'u006_003'; %B5 tperiod
% trxexp = 'u006_008'; %B2 tperiod

% anim = 'rm1';
% tempexp = 'u003_001'; %Experiment that gets aligned "to"
% % trxexp = 'u003_003'; %B5 tperiod theta [0 90]; Experiment that gets transformed
% % trxexp = 'u003_004'; %B5 tperiod theta [45 135]
% % trxexp = 'u003_006'; %B2 tperiod theta [0 90]
% % trxexp = 'u003_007'; %B2 tperiod theta [45 135]

% anim = 'rl7';
% tempexp = 'u006_000'; %Experiment that gets aligned "to"
% % trxexp = 'u006_002'; %B5 tperiod theta [0 90]; Experiment that gets transformed
% % trxexp = 'u006_003'; %B5 tperiod theta [45 135]
% % trxexp = 'u006_005'; %B1 tperiod theta [0 90]
% % trxexp = 'u006_006'; %B1 tperiod theta [45 135]

% anim = 'rm5';
% tempexp = 'u003_000'; %Experiment that gets aligned "to"
% % trxexp = 'u003_003'; %B1 sfreq theta [0 90]; Experiment that gets transformed
% % trxexp = 'u003_004'; %B1 sfreq  theta [45 135]
% % trxexp = 'u003_007'; %B5 sfreq  theta [0 90]
% trxexp = 'u003_008'; %B5 sfreq  theta [45 135]

% anim = 'rm4';
% tempexp = 'u003_000'; %Experiment that gets aligned "to"
% % trxexp = 'u003_003'; %B1 sfreq theta [0 90]; Experiment that gets transformed
% % trxexp = 'u003_004'; %B1 sfreq  theta [45 135]
% % trxexp = 'u003_006'; %B5 sfreq  theta [0 90]
% % trxexp = 'u003_007'; %B5 sfreq  theta [45 135]

% anim = 'rm0';
% tempexp = 'u003_007'; %Experiment that gets aligned "to"
% trxexp = 'u003_002'; %B1 sfreq theta [0 90]; Experiment that gets transformed
% % trxexp = 'u003_003'; %B1 sfreq  theta [45 135]
% % trxexp = 'u003_005'; %B5 sfreq  theta [0 90]
% % trxexp = 'u003_006'; %B5 sfreq  theta [45 135]

% anim = 'rw3';
% tempexp = 'u001_001'; %Experiment that gets aligned "to"
% % trxexp = 'u001_005'; %B5 tperiod theta [0 90]; Experiment that gets transformed
% % trxexp = 'u001_006'; %B5 tperiod  theta [45 135]
% % trxexp = 'u001_008'; %B5 atropine tperiod  theta [0 90]
% % trxexp = 'u001_013'; %B5 atropine tperiod  theta [45 135]

% anim = 'rs9';
% tempexp = 'u002_010'; %Experiment that gets aligned "to"
% % trxexp = 'u002_003'; %B5 tperiod theta [0 90]; Experiment that gets transformed
% % trxexp = 'u002_004'; %B5 tperiod  theta [45 135]
% % trxexp = 'u002_007'; %B5 atropine tperiod  theta [0 90]
% % trxexp = 'u002_009'; %B5 atropine tperiod  theta [45 135]

% anim = 'rr1';
% tempexp = 'u003_008'; %Experiment that gets aligned "to"
% % trxexp = 'u003_002'; %B5 tperiod theta [0 90]; Experiment that gets transformed
% % trxexp = 'u003_003'; %B5 tperiod  theta [45 135]
% % trxexp = 'u003_006'; %B5 atropine tperiod  theta [0 90]
% % trxexp = 'u003_007'; %B5 atropine tperiod  theta [45 135]

% anim = 'rr0';
% tempexp = 'u003_012'; %Experiment that gets aligned "to"
% % trxexp = 'u003_006'; %B5 tperiod theta [0 90]; Experiment that gets transformed
% % trxexp = 'u003_007'; %B5 tperiod  theta [45 135]
% % trxexp = 'u003_009'; %B5 atropine tperiod  theta [0 90]
% % trxexp = 'u003_010'; %B5 atropine tperiod  theta [45 135]

% anim = 'rm4';
% tempexp = 'u005_000'; %Experiment that gets aligned "to"
% % trxexp = 'u005_002'; %B5 tperiod theta [0 90]; Experiment that gets transformed
% % trxexp = 'u005_003'; %B5 tperiod  theta [45 135]
% % trxexp = 'u005_005'; %B5 atropine tperiod  theta [0 90]
% % trxexp = 'u005_006'; %B5 atropine tperiod  theta [45 135]

anim = 'se3';
tempexp = 'u001_002'; %Experiment that gets aligned "to"
% trxexp = 'u001_004'; %B5 sfreq theta [0 90]; Experiment that gets transformed
trxexp = 'u001_005'; %B5 sfreq  theta [45 135]

% anim = 'se4';
% tempexp = 'u001_001'; %Experiment that gets aligned "to"
% % trxexp = 'u001_003'; %B5 sfreq theta [0 90]; Experiment that gets transformed
% % trxexp = 'u001_004'; %B5 sfreq  theta [45 135]


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



imloc = [dataroot anim '\' tempexp '\' tempexp '_000_f50.mat']; %pick a random frame
S = load(imloc);  %Returns the contents in the .mat under the structure S
im_output = double(S.im);

im_output = im_output-prctile(im_output(:),1);
im_output(find(im_output<0)) = 0;
im_output = im_output/prctile(im_output(:),99.9);
im_output(find(im_output>1)) = 1;

disp(['Aligning to template'])

imloc = [dataroot anim '\' trxexp '\' trxexp '_000_f50.mat'];
S = load(imloc);  %Returns the contents in the .mat under the structure S
im_input = double(S.im);

im_input = im_input-prctile(im_input(:),1);
im_input(find(im_input<0)) = 0;
im_input = im_input/prctile(im_input(:),99.9);
im_input(find(im_input>1)) = 1;

[movingPoints,fixedPoints] = cpselect(im_input,im_output,'Wait',true);

trx = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');

im_in = im_input;

%%
path = [spath 'align_' anim '_' trxexp(2:end) '_to_' tempexp(2:end)];
uisave('trx',path)


%%

[xp yp] = meshgrid(1:20:size(im_output,2),1:20:size(im_output,1));
figure
subplot(1,2,1)
imagesc(im_output), colormap gray,axis image
hold on, plot(xp(:),yp(:),'.r')
title('template')

Ref_Im = imref2d(size(im_output));
dum = imwarp(im_in,trx,'OutputView',Ref_Im);
subplot(1,2,2)
imagesc(dum), colormap gray,  axis image

hold on, plot(xp(:),yp(:),'.r')

