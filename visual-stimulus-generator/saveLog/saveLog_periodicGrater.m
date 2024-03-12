function saveLog_periodicGrater


global Mstate loopTrial TimingInfo

%root = '/Matlab_code/log_files/';

root = '~/Desktop/log_files/';

%rootnet = ['/Volumes/neurostuff/log_files/' Mstate.anim '/'];

expt = [Mstate.anim '_' Mstate.unit '_' Mstate.expt];

fname = [root expt '.mat'];
%fnamenet = [rootnet expt '.mat'];

frate = Mstate.refresh_rate;

%%%

if loopTrial == 1
    save(fname,'frate')
    %save(fnamenet,'frate')    
end

save(fname,'TimingInfo','-append')  %This will just write over the previous version on each trial
%save(fnamenet,'TimingInfo','-append')


