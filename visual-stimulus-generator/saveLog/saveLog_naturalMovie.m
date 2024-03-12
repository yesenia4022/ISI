function saveLog_naturalMovie(seqs)


global Mstate loopTrial TimingInfo imageOnset


%root = '/Matlab_code/log_files/';

root = '~/Desktop/log_files/';

%rootnet = ['/Volumes/neurostuff/log_files/' Mstate.anim '/'];

expt = [Mstate.anim '_' Mstate.unit '_' Mstate.expt];

fname = [root expt '.mat'];
%fnamenet = [rootnet expt '.mat'];

frate = Mstate.refresh_rate;

%%%

if isempty(seqs) %Indicates that it is called from 'play'
    % randomized image ID's
    randImageID = ['randlog_T' num2str(loopTrial)];
    eval([randImageID '.seqs = seqs;'])
    
    % time stamps for each image onset
    stimOnset  = ['stimOnset_', num2str(loopTrial)];
    eval([stimOnset '.timestamps = imageOnset;']);
%     
    %save(fname,'TimingInfo','imageOnset');  %This will just write over the previous version on each trial
    %save(fnamenet,'TimingInfo','-append')

     if loopTrial == 1
%         save(fname,['randlog_T' num2str(loopTrial)])
%         save(fname,'frate',['stimOnset_',num2str(loopTrial)], '-append')
        save(fname, 'frate',['randlog_T' num2str(loopTrial)],['stimOnset_',num2str(loopTrial)]);
    else
        save(fname,['randlog_T' num2str(loopTrial)],['stimOnset_',num2str(loopTrial)],'-append')
%         save(fname, ['stimOnset_',num2str(loopTrial)], '-append');
     end
    
     
else %called from make
    
    % randomized image ID's
    randImageID = ['randlog_T' num2str(loopTrial)];
    eval([randImageID '.seqs = seqs;'])
    
    % time stamps for each image onset
    stimOnset  = ['stimOnset_', num2str(loopTrial)];
    eval([stimOnset '.timestamps = imageOnset;']) ;
    
    if loopTrial == 1
%         save(fname,['randlog_T' num2str(loopTrial)])
%         save(fname,'frate',['stimOnset_',num2str(loopTrial)], '-append')
        save(fname, 'frate',['randlog_T' num2str(loopTrial)], ['stimOnset_',num2str(loopTrial)]);
    else
        save(fname,['randlog_T' num2str(loopTrial)],['stimOnset_',num2str(loopTrial)],'-append')
%         save(fname, ['stimOnset_',num2str(loopTrial)], '-append');
    end
    
end


