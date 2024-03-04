function ManualCellSorter(im)

global posall

posall = [0 0];

figure, imagesc(im), colormap gray

axis image;

fh = gcf;

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

global posall

pos = round(get(event_obj,'Position')); %pos(1) is column dimension

hold on, plot(pos(1),pos(2),'.r')

d = sqrt((pos(1)-posall(end,1))^2 + (pos(2)-posall(end,2))^2);
if d > 3
    posall = [posall; pos];
end

