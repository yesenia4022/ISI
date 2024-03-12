function altazimuth = altazimuthSwitch(ori)

%N.B. the screen is generally rotated for the mouse
% if ori == 0 || ori == 180    
%     altazimuth = 'altitude';    
% elseif ori == 90 || ori == 270    
%     altazimuth = 'azimuth';    
% end

% % When screen horizontal
% if ori == 0 || ori == 180    
%     altazimuth = 'azimuth';    
% elseif ori == 90 || ori == 270    
%     altazimuth = 'altitude';    
% end
% When screen horizontal & I'm trying different ori combos: GCR 6/01/20
if ori == 0 || ori == 180    
    altazimuth = 'azimuth';    
elseif ori ~= 0 || ori ~= 180
    altazimuth = 'altitude';    
end


