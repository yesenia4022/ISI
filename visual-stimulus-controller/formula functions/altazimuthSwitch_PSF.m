function [AZM1 AZM2] = altazimuthSwitch_PSF(ori,ori2)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ori == 0 | ori == 180
%     AZM1 = 'azimuth'
% elseif ori == 90 | ori == 270
%     AZM1 = 'altitude';
% end
% 
% if ori2 == 0 | ori2 == 180
%     AZM2 = 'azimuth';
% elseif ori2 == 90 | ori2 == 270
%     AZM2 = 'altitude';
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gaby edits for trying Declan's non stand ori
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6/01/20
if ori == 0 | ori == 180
    AZM1 = 'altitude' %8/16/21 edit to match Kaltsky
%   AZM1 = 'azimuth'

elseif ori ~= 0 | ori ~= 180
    AZM1 = 'altitude';
end

if ori2 == 0 | ori2 == 180
    AZM2 = 'azimuth';
elseif ori2 == 90 | ori2 == 270
    AZM2 = 'altitude';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%