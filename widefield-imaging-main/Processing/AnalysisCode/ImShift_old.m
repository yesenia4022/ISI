function imout = ImShift(im,mbest,nbest)

%'mbest' and 'nbest' is the amount to shift 'im' in the first and second dimension,
%respectively.
%mbest and nbest are acquired with 'getShiftVals.m'

dim = size(im);

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
