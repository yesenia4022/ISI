function saveLog_rain(domains,seq)

global Mstate loopTrial TimingInfo

%root = '/Matlab_code/log_files/';

root = '~/Desktop/log_files/';
%rootnet = ['/Volumes/neurostuff/log_files/' Mstate.anim '/'];

        
expt = [Mstate.anim '_' Mstate.unit '_' Mstate.expt];

fname = [root expt '.mat'];
%fnamenet = [rootnet expt '.mat'];

frate = Mstate.refresh_rate;

if isempty(seq) %Indicates that it is called from 'play'
    
    save(fname,'TimingInfo','-append')  %This will just write over the previous version on each trial
    %save(fnamenet,'TimingInfo','-append')
    
else  %from 'make'
    
    if loopTrial == 1  %save domains only on first trial
        
        if ~exist(fname)
            save(fname,'domains','frate')
            %save(fnamenet,'domains','frate')
        end
        
    end
    
    %Append file with sequence info on every trial
    
    eval(['rseed' num2str(loopTrial) '=seq;' ])
    eval(['save ' fname ' rseed' num2str(loopTrial) ' -append'])
    %eval(['save ' fnamenet ' rseed' num2str(loopTrial) ' -append'])
    
end
