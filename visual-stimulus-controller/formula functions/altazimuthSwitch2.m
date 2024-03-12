function altazimuth2 = altazimuthSwitch2(ori2)

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
if ori2 == 0 || ori2 == 180    
    altazimuth2 = 'azimuth';    
elseif ori2 ~= 0 || ori2 ~= 180
    altazimuth2 = 'altitude';    
end


