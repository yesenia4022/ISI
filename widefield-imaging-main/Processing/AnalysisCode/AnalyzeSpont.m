function AnalyzeSpont

global CHs bwCell1 bwCell2

filepath = 'c:\2p_data\gm0\Spont0.tif';
%filepath = 'c:\2p_data\gm0\GFPpic005.tif';

tf = imformats('tif');
info = feval(tf.info, filepath);

Nimages = length(info);

tf = imformats('tif');

Nchan = 2;
shiftflag = -1;

%%%%%Get Data%%%%%
framestart = 1;
framestop = Nimages;
clear CHs
chvec = [1 1 1];
i = 1;
if chvec(1) == 1
    k = 1;
    for frame=framestart:Nchan:framestop
        A = feval(tf.read, filepath, frame);
        CHs{i}(:,:,k) = double(A);
        k = k+1;
    end
    i = 2;
end

if chvec(2) == 1
    k = 1;
    for frame=(framestart+1):Nchan:framestop
        A = feval(tf.read, filepath, frame);
        CHs{i}(:,:,k) = double(A);
        k = k+1;
    end
    i = i+1;
end
%%%%%%%%%%%%%%%

dim = size(CHs{1});

%%%Apply Shift correction

if shiftflag
    temp1 = mean(CHs{1},3);
    temp2 = mean(CHs{2},3);  %set template image
    for z = 1:length(CHs{1}(1,1,:))
        imdum = CHs{1}(:,:,z);
        [mbest nbest] = getShiftVals(imdum,temp1);  %get the transformation
        CHs{1}(:,:,z) = circshift(CHs{1}(:,:,z),[-mbest -nbest]); %transform
        
        imdum = CHs{2}(:,:,z);
        [mbest nbest] = getShiftVals(imdum,temp2);  %get the transformation
        CHs{2}(:,:,z) = circshift(CHs{2}(:,:,z),[-mbest -nbest]); %transform
    end
end

% temp = (CHs{1}(:,:,2) + CHs{2}(:,:,2))/2;  %set template image
% 
% for z = 1:length(CHs{1}(1,1,:))
%     imdum = (CHs{1}(:,:,z) + CHs{2}(:,:,z))/2;
%     [mbest nbest] = getShiftVals(imdum,temp);  %get the transformation
% 
%     CHs{1}(:,:,z) = ImShift(CHs{1}(:,:,z),mbest,nbest);  %transform
%     CHs{2}(:,:,z) = ImShift(CHs{2}(:,:,z),mbest,nbest);
%     
%     dum = CHs{1}(:,:,z);
%     normer = prctile(dum(:),10);
%     CHs{1}(:,:,z) = (CHs{1}(:,:,z)-normer)/normer;  
%     dum = CHs{2}(:,:,z);
%     normer = prctile(dum(:),10);
%     CHs{2}(:,:,z) = (CHs{2}(:,:,z)-normer)/normer;  
%     z
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Make the mask%%%%%%%%%%%
CH1 = 0; CH2 = 0;

temp1 = sum(CHs{1},3);
temp2 = sum(CHs{2},3);

[dum bwCell1] = LocalZ(temp1,14,.9);

mask = [0 1 0;1 1 1; 0 1 0];
bwCell1 = imopen(bwCell1,mask);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear xloc yloc cellTcourse1 cellTcourse2
cellid = bwlabel(bwCell1);
uid = unique(cellid);

cellTcourse = zeros(length(uid),length(CHs{1}(1,1,:)));

for z = 1:dim(3)
    imdum = CHs{1}(:,:,z);
    for i = 1:length(uid)
        idx = find(cellid(:) == uid(i));
        xlocdum = ceil(idx/dim(1));
        ylocdum = idx-(xlocdum-1)*dim(1);
        xloc(i) = round(mean(xlocdum));
        yloc(i) = round(mean(ylocdum));
        
        
        cellTcourse(i,z) = mean(imdum(idx));
    end
end

cellTcourse = zscore(cellTcourse')';
covMat = (cellTcourse*cellTcourse')/length(cellTcourse(1,:));

figure,imagesc(covMat), colorbar

clear cc D
k = 1;
for i = 1:length(uid)
    for j = 1:length(uid)
        Dx = abs(xloc(i)-xloc(j));
        Dy = abs(yloc(i)-yloc(j));

        D(i,j) = sqrt((xloc(i)-xloc(j))^2 + (yloc(i)-yloc(j))^2);

        cc(i,j) = covMat(j,i);
        k = k+1;
    end
end

clear R
for i = 1:length(uid)
    Ddum = D(:,i);
    id = find(Ddum ~= 0 & ~isnan(Ddum));
    Rdum = corrcoef(D(id,i),covMat(id,i))
    R(i) = Rdum(2,1);
end

[dum i] = max(R);
Ddum = D(:,i);
id = find(Ddum ~= 0);
figure,scatter(D(id,i)*400/256,covMat(id,i),'.')

im = zeros(dim(1),dim(2));
for i = 1:length(uid)
    idx = find(cellid(:) == uid(i));
    im(idx) = R(i);
end
figure,imagesc(im)

figure,scatter(D(:),cc(:),'.')
corrcoef(D(:),cc(:));

%%%%%

Nc = length(uid);
id = find(cellid ~= 0);
mask = zeros(dim(1),dim(2));
mask(id) = 1;

figure
for i = 1:Nc
    im = zeros(dim(1),dim(2));
    for j = 1:Nc
        idx = find(cellid(:) == uid(j));
        im(idx) = cc(i,j);
    end
    subplot(floor(sqrt(Nc)),ceil(sqrt(Nc)),i)
    imagesc(im,'AlphaData',mask,[-1 1]), axis off
end



%Cellurize(sfmap1,bwCell1);
