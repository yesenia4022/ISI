function [mbest nbest] = getShiftVals(im,temp)

% dim = size(Tens);
% hh = hann(dim(1))*hann(dim(2))';
% hh = hh/sum(hh);
% delta = zeros(dim(1),dim(2));
% [ma id] = max(hh(:));
% delta(id(1)) = 1;
% hh = delta-hh; %high pass

% for i = 1:dim(3)
%     dum = squeeze(Tens(:,:,i));
%     %dum = ifft2(abs(fft2(hh)).*fft2(dum));
%     Tens(:,:,i) = (dum-mean(dum(:)))/std(dum(:));
% end

W = 3;
shiftdom = -W:W;
dim = size(im);

y = 1;
for m = shiftdom(1):shiftdom(end)
    z = 1;
    for n = shiftdom(1):shiftdom(end)

        [m1 m2 n1 n2] = getrangePos(m,n,dim(1),dim(2));
        temppiece = temp(m1:m2,n1:n2);

        [m1 m2 n1 n2] = getrangeNeg(m,n,dim(1),dim(2));
        impiece = im(m1:m2,n1:n2);

        R = corrcoef(temppiece(:),impiece(:));
        cc(y,z) = R(1,2);
        z = z+1;
    end
    y = y+1;
end

[idr idc] = find(cc == max(cc(:)));
mbest = shiftdom(idr);
nbest = shiftdom(idc);

[m1 m2 n1 n2] = getrangeNeg(mbest,nbest,dim(1),dim(2));
impiece = im(m1:m2,n1:n2);

imout = NaN*zeros(size(im));
[m1 m2 n1 n2] = getrangePos(mbest,nbest,dim(1),dim(2));
imout(m1:m2,n1:n2) = impiece;
    



function [m1 m2 n1 n2] = getrangePos(m,n,rows,cols)
%for the stationary object, e.g. the template

m1 = m+1;
n1 = n+1;
m1 = max(1,m1);
n1 = max(1,n1);

m2 = rows+m;
n2 = cols+n;
m2 = min(rows,m2);
n2 = min(cols,n2);


function [m1 m2 n1 n2] = getrangeNeg(m,n,rows,cols)
%for the shifting object, e.g. the images getting registered to the
%template

m1 = -m+1;
n1 = -n+1;
m1 = max(1,m1);
n1 = max(1,n1);

m2 = rows-m;
n2 = cols-n;
m2 = min(rows,m2);
n2 = min(cols,n2);
