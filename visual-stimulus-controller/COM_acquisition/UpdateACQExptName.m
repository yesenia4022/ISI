function UpdateACQExptName

global Mstate TwoPcomState

% set(0,'DefaultTextFontName','helvetica','DefaultTextFontAngle','normal','DefaultTextColor',[0 0 0])
% button = questdlg(sprintf('Are you sure you want to save the data\nand advance to the next experiment?'));
% set(0,'DefaultTextFontName','helvetica','DefaultTextFontAngle','oblique','DefaultTextColor',[1 1 0])

if Mstate.twoP
    
    %Send to Scanbox
    fprintf(TwoPcomState.serialPortHandle,['A' Mstate.anim]);
    fprintf(TwoPcomState.serialPortHandle,['U' Mstate.unit]);
    fprintf(TwoPcomState.serialPortHandle,['E' Mstate.expt]);
    
end


%if Mstate.WF  %I don't like this old conditional because somebody could
%increment the experiment while its unclicked, click it, and then run, thus
%leaving a mismatch. The conditional below matches them up if imager is
%simply open.
if ~isempty(findobj('Tag','animaltxt')) %This is a field of the 'imager'
    
    %Update the 'imager' fields, even though they are not used in identifying
    %file save location (the Stimulator fields are)
    set(findobj('Tag','animaltxt'),'String',Mstate.anim);
    set(findobj('Tag','expttxt'),'String',Mstate.expt);
    set(findobj('Tag','unittxt'),'String',Mstate.unit);
    
    set(findobj('Tag','ISIdataRoot'),'String',Mstate.dataRoot);
    set(findobj('Tag','analyzerRoot'),'String',Mstate.analyzerRoot);  
    
end



if Mstate.Ephys
    
    
end