function updateTempReading

global imager_handle ACQ

%src = getselectedsource(ACQ.vid);
%temp = round(src.DeviceTemperature*100)/100;
%set(imager_handle.cameraTemp,'String',['Camera Temp: ' num2str(temp) ' deg']);