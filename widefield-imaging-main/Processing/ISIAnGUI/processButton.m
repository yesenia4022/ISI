function processButton

global f0m f0m_var Tens Tens_var AUE bsflag Flim repDom G_handles

bsflag = 0;

t0 = cputime;

bsflag = get(G_handles.basesub,'Value');
Flim = str2double(get(G_handles.epistart,'String'));  %Frame start in ms (to average)
Flim(2) = str2double(get(G_handles.epistop,'String')); %Frame stop in ms (to average)
b = str2double(get(G_handles.bstart,'String')); %in msec as well
b(2) = str2double(get(G_handles.bstop,'String'));

slowMo = get(G_handles.slowMotionFlag,'Value');  %lateral movement correction

repdum = get(G_handles.repDom,'string');
if strcmp(repdum,'All')
    repDom = 1:getnorepeats(1);
else
    eval(['repDom = [' repdum '];'])
end

setacqinfo(1)  %need to set to first trial so that it doesn't use a blank trial when doing CondTensor

set(G_handles.status,'string','Processing...'), drawnow

if get(G_handles.F0flag,'value') %Get F0 images (and possibly cell mask traces)
    [Tens Tens_var] = CondTensor3(b,slowMo);  %%Compute entire space time block for each condition
    varflag = 1;
    
    
    if ~get(G_handles.F1flag,'Value')
        f0m = CondF0(Tens,Flim); %%%Compute |F1| over time interval [Flim(1) Flim(2)]%%%
        if ~isempty(Tens_var{1})
            f0m_var = CondF0(Tens_var,Flim); %%%Compute |F1| over time interval [Flim(1) Flim(2)]%%%
        end
    else
        f0m = CondF1(Tens,Flim); %%%Compute mean over time interval [Flim(1) Flim(2)]%%%
        if ~isempty(Tens_var{1})
            f0m_var = CondF1(Tens_var,Flim); %%%Compute |F1| over time interval [Flim(1) Flim(2)]%%%
        end
    end
    
    %f0m_var = CondF0(Tens_var,Flim);
else %only get cell mask traces (for all trials)
    CondMaskdata(slowMo,fastMo)
end


set(G_handles.status,'string','Done'), drawnow
%sound(.6*sin(2*pi*400/(1200)*(0:400)),1200)  %Signal done

t1 = cputime-t0;

set(G_handles.time,'string',num2str(t1))

set(G_handles.loaded,'string',AUE)

set(G_handles.plot,'enable','on')

set(G_handles.save,'enable','on')
