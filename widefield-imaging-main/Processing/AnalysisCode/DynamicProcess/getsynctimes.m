function synctimes = getsynctimes(x)

global ACQinfo

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

mid = (max(x)+min(x))/2;
x = sign(x'-mid)/2;
%x = diff(-x);
x = -x;
x = x(2:end)-x(1:end-1);
id = find(x == 1);
synctimes = id*ptime;  %There will be jitter error around flyback time, at most (i.e. <.5ms)

dsync = diff(synctimes);
flashperiod = min(dsync);
id = find(dsync>flashperiod*1.5);

while ~isempty(id)
    
    synctimes = [synctimes(1:id(1)) synctimes(id(1))+flashperiod synctimes(id(1)+1:end)];
    
    dsync = diff(synctimes);
    id = find(dsync>flashperiod*1.5);
    
end
