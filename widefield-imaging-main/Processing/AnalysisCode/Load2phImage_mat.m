function [CHs] = Load2phImage_mat(frame,chvec,varargin)

%varargin can be the trial number (1 to N)
%if you want to use the preset condition/repeat to define the trial, leave out varargin.

global twophDATADIR AUE


trial = varargin{1};

filepath = [twophDATADIR AUE ' ' sprintf('%03d',trial) '.mat'];

%load(filepath,'state')
load(filepath)

if chvec(1) == 1    
    CHs{1} = double(acquiredData{1}(:,:,frame));    
end

if chvec(2) == 1    
    CHs{2} = double(acquiredData{2}(:,:,frame));     
end

if chvec(3) == 1    
    CHs{3} = double(acquiredData{3}(:,:,frame));     
end

if chvec(4) == 1    
    CHs{4} = double(acquiredData{4}(:,:,frame));     
end


% load(filepath,'acquiredData')
% 
% if chvec(1) == 1    
%     CHs{1} = double(acquiredData{1}(:,:,frame));    
% end
% 
% if chvec(2) == 1    
%     CHs{2} = double(acquiredData{2}(:,:,frame));     
% end
% 
% if chvec(3) == 1    
%     CHs{3} = double(acquiredData{3}(:,:,frame));     
% end
% 
% if chvec(4) == 1    
%     CHs{4} = double(acquiredData{4}(:,:,frame));     
% end

