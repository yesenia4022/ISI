function [imMat gdom bdom] = getBGmatrix

global f0m Analyzer

odom = getdomain('ori');
gdom = getdomain('greengain');
bdom = getdomain('bluegain');

imMat = cell(length(bdom),length(gdom));
for cond = 1:length(f0m)-1
    
    greengain = getVal('greengain',cond);
    gid = find(gdom == greengain);
    
    bluegain = getVal('bluegain',cond);
    bid = find(bdom == bluegain);
    
    if isempty(imMat{bid,gid})
        imMat{bid,gid} = 0;
    end
    
    imMat{bid,gid} = imMat{bid,gid}+f0m{cond}/length(odom);
    
    ggain(bid,gid) = greengain;
    bgain(bid,gid) = bluegain;
    
end

hh = fspecial('gaussian',size(imMat{1,1}),2);
imMatF = cell(1,size(imMat,1));
for i = 1:size(imMat,1)
    for j = 1:size(imMat,2)        
        imMatF{i}(:,:,j) = ifft2(abs(fft2(hh)).*fft2(imMat{i,j}));        
    end
end

figure
for i = 1:length(bdom)
    dum = -imMatF{i} + repmat(max(imMatF{i},[],3),[1 1 length(gdom)]);
    dum = dum./repmat(sum(dum,3),[1,1,length(gdom)]);
    
    mag = mean(imMatF{i},3);
    mi = prctile(mag(:),1);    
    mag = mag-mi;
    mag(find(mag<0)) = 0;
    
    ma = prctile(mag(:),50);
    mag = mag/ma;
    mag(find(mag>1)) = 1;  
    
    
    CM = 0;
    for j = 1:length(gdom)
        
        CM = CM + gdom(j)*dum(:,:,j);
        
    end
    
    mi = prctile(CM(:),3);
    ma = prctile(CM(:),97);
    
    subplot(1,3,i),imagesc(CM,'AlphaData',mag,[-.6 .6]), axis image, colorbar
    colormap jet
    
end

Miso = (imMatF{1}(:,:,6)+imMatF{1}(:,:,2))/2;
Uiso = imMatF{2}(:,:,2);
figure,imagesc((Miso-Uiso)./(Miso+Uiso),[-.8 .8])
    

function val = getVal(pstring,cond)

global Analyzer

Nparam = length(Analyzer.loops.conds{1}.val);

for k = 1:Nparam
    if strcmp(Analyzer.loops.conds{cond}.symbol{k},pstring)
        val = Analyzer.loops.conds{cond}.val{k};
        break
    end
end

