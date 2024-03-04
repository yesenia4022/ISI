function SaveTrial(trialno)

global Tens fname ACQ

trial = trialno-1;  %trial number in file name starts with n = 0

%These used to look to the imager GUI fields.  
%For safety and general ease of keeping things straight, I now have it
%use the fields of "Stimulator"
animal = get(findobj('Tag','animal'),'String');
unit   = get(findobj('Tag','unitcb'),'String');
expt   = get(findobj('Tag','exptcb'),'String');
datadir= get(findobj('Tag','dataRoot'),'String');

dd = [datadir '\' lower(animal) '\u' unit '_' expt];

fname = sprintf('%s\\u%s_%s',dd,unit,expt);
fname = [fname  '_' sprintf('%03d',trial)];

%Enter saving loop

tsize = size(Tens,3);

% if ~isempty(ACQ.ROIcrop)
%     xmin = ACQ.ROIcrop(1);
%     xmax = ACQ.ROIcrop(1)+ACQ.ROIcrop(3)-1;
%     ymin = ACQ.ROIcrop(2);
%     ymax = ACQ.ROIcrop(2)+ACQ.ROIcrop(4)-1;
% else   %Default is for this to be empty and save the entire image
%     xmin = 1;
%     xmax = size(Tens,2);
%     ymin = 1;
%     ymax = size(Tens,1);
% end

xmin = 1;
xmax = size(Tens,2);
ymin = 1;
ymax = size(Tens,1);

% delete this after 6/26
ACQ.offLineSpatialBinning = 1;
D = ACQ.offLineSpatialBinning;
xmin = ceil(xmin/D);
ymin = ceil(ymin/D);
xmax = ceil(xmax/D);
ymax = ceil(ymax/D);

for n = 1:tsize
    im = Tens(ymin:ymax,xmin:xmax,n);
    var = ['f' num2str(n)];
    fnamedum = [fname '_' var];
    save(fnamedum,'im')
end


