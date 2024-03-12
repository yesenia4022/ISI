function saveLog_uncertaintyStimulus(beforeOrAfter)


global Mstate loopTrial 

global stimulusFolder scrambledFileNames scrambledContrasts stimOnset stimOffset 

%root = '/Matlab_code/log_files/';

root = '~/Desktop/log_files/';

%rootnet = ['/Volumes/neurostuff/log_files/' Mstate.anim '/'];

expt = [Mstate.anim '_' Mstate.unit '_' Mstate.expt];

fname = [root expt '.mat'];
%fnamenet = [rootnet expt '.mat'];

%frate = Mstate.refresh_rate;
BEFORE  = 0;
AFTER   = 1;
if(beforeOrAfter == BEFORE)
    save(fname, 'stimulusFolder', ...
    'scrambledFileNames', 'scrambledContrasts');
elseif(beforeOrAfter == AFTER)
    if(exist(fname,'file') > 0)
        save(fname, 'stimOnset', 'stimOffset', '-append');
    else
        save(fname, 'stimulusFolder', ...
        'scrambledFileNames', 'scrambledContrasts', ...
        'stimOnset', 'stimOffset');
    end
else
    save(fname, 'stimulusFolder', ...
    'scrambledFileNames', 'scrambledContrasts', ...
    'stimOnset', 'stimOffset');
end
%{
if isempty(seqs) %Indicates that it is called from 'play'
    
    save(fname,'TimingInfo','-append')  %This will just write over the previous version on each trial
    %save(fnamenet,'TimingInfo','-append')

else %called from make
    
    basename = ['randlog_T' num2str(loopTrial)];
    
    eval([basename '.seqs = seqs;'])
    
    if loopTrial == 1
        save(fname,['randlog_T' num2str(loopTrial)])
        %save(fnamenet,['randlog_T' num2str(loopTrial)])
        
        save(fname,'frate','-append')
        %save(fnamenet,'frate','-append')
        
    else
        save(fname,['randlog_T' num2str(loopTrial)],'-append')
        %save(fnamenet,['randlog_T' num2str(loopTrial)],'-append')
        
    end
    
end

%}
