function [flag dd] = checkforOverwrite

global running

%Use the fields of the Stimulator, not the 'imager'. The code looks to the Stimulator fields
%to actually id the save location.  Although, they (imager/Stimulator fields) should be the same, of course.
%Its just for safety.
animal = get(findobj('Tag','animal'),'String');
unit   = get(findobj('Tag','unitcb'),'String');
expt   = get(findobj('Tag','exptcb'),'String');
datadir= get(findobj('Tag','dataRoot'),'String');

dd = [datadir '\' animal '\u' unit '_' expt]

flag = 0;

if(exist(dd))
    warndlg('Directory exists!!!  Make sure you are not overwritting old data!  Please check current animal, unit and expt values.  I will now abort this sampling request.','!!! Warning !!!')
    flag = 1;
end
