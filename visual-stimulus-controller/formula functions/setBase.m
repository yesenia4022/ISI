function gb = setBase(base_scale,gb)
% 'base_scale' inputs for each experiment: [.03125 .0625 .125 .25 .5]
gb = base_scale*2*gb;

gb(3) = base_scale; % RGB bases
gb(4) = base_scale*2*128; % Background

% In Looper window, formula:
% gb = getGreenBlueGain(theta,1,1); gb = setBase(0.5,gb); greengain = gb(1); bluegain = gb(2); redbase = gb(3); greenbase = gb(3); bluebase = gb(3); background = gb(4); 

end