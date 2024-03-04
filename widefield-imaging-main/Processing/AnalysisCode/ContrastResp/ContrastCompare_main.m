% load expt
% process raw



global f0m maskS

%set roi

%bw = ones(256,256);

%bwCell1 = MakeCellMask(14,.8,3);

hh = [];
[f0 ori cont] = GetOrivsContrast(f0m,hh);

PopContrast(f0,ori,cont,maskS.bwCell1);





%%%%%%
%%%%%%
%fc4
load('F:\neurostuff\Masks\fc4exp1')

