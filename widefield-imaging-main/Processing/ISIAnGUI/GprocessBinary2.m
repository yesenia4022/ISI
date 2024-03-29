function [binarymap binmask] = GprocessBinary2(f0,hh,hnorm,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2 is used for the SF/OD paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global bsflag Analyzer symbolInfo G_handles
%Each element of the cell array 'f0dum' is the average image for the
%corresponding condition

idsym = symbolInfo.ID(1);

idsym2 = symbolInfo.ID(2);

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
if bflag
    f0blank = f0{end};
    f0(end) = [];
end

for i = 1:length(f0)
    bindomcond(i) = Analyzer.loops.conds{i}.val{idsym};
    
    secdomcond(i) = Analyzer.loops.conds{i}.val{idsym2};
end
bindom = unique(bindomcond);
secdom = unique(secdomcond);  %domain of secondary parameter... only useful if ~empty(varargin)

dim = size(f0{1});
Tens = zeros(dim(1),dim(2),length(f0),'single'); %preallocate
for k = 1:length(f0)
    Tens(:,:,k) = f0{k};
end

%I used to normalize by the blank if baseline subtraction was not checked.
%This is now done on the Tensor beforehand in CondTensor3

%If I subtracted the initial baseline, it doesn't help to normalize by blank, it just adds noise...
% mi = min(Tens,[],3);
% for k = 1:length(Tens(1,1,:))
%     if bflag == 1  %if a blank exists
%         Tens(:,:,k) = Tens(:,:,k)-f0blank;
%     else
%         Tens(:,:,k) = Tens(:,:,k)-mi;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(varargin)  %limit to a subsection of values within the second parameter domain
    secdomLimit = varargin{1};
    im1 = 0; counter = 0;
    for j = 1:length(secdomLimit)
        id = find(secdomcond == secdomLimit(j) & bindomcond == bindom(1));
        im1 = im1 + nansum(Tens(:,:,id),3);
        counter = counter + length(id);
    end
    im1 = im1/counter;
    
    im2 = 0; counter = 0;
    for j = 1:length(secdomLimit)
        id = find(secdomcond == secdomLimit(j) & bindomcond == bindom(2));
        im2 = im2 + sum(Tens(:,:,id),3);
        counter = counter + length(id);
    end
    im2 = im2/counter;

else  %ignore the other domains and average over everything

    id = find(bindomcond == bindom(1));
    im1 = mean(Tens(:,:,id),3);
    id = find(bindomcond == bindom(2));
    im2 = mean(Tens(:,:,id),3);

end

idN = find(isnan(im1.*im2));
im1(idN) = 0; im2(idN) = 0;

im1 = phi(im1); im2 = phi(im2);
binarymap = im2-im1;   %keep it as 2-1, to be consistent with eye logic
normer = im1+im2;


%if a filter exists, use it...
id = find(isnan(binarymap));
binarymap(id) = 0;
normer(id) = 0;
if ~isempty(hh)
    binarymap = ifft2(abs(fft2(hh)).*fft2(binarymap));
    normer = ifft2(abs(fft2(hnorm)).*fft2(normer));
end

%Normalize after smoothing to make division more stable
binarymap = binarymap./normer;

dum = im1+im2;
dum = medfilt2(dum,[5 5]);
binmask = sign(dum);
se = strel('disk',1);
binmask = imerode(binmask,se);

set(G_handles.Lwidth,'string','1');
h = makeMapFilter;
binmask = ifft2(fft2(double(binmask)).*abs(fft2(h)));