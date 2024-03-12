function prepACQ

%Prep Scanimage... give it trial time and hit the loop button

total_time = getParamVal('stim_time',1);
predelay = getParamVal('predelay',1);
postdelay = getParamVal('postdelay',1);

%Send trial length values
Stimulus_localCallback(['StimTimeInfo=' 13 ...
    num2str(total_time) 13 ...
    num2str(predelay) 13 ...
    num2str(postdelay)]);


%Tell ScanImage to hit Loop button and then wait for response:
Stimulus_localCallback('loop'); %Tell ScanImage to start 'loop' button

