function configPstate_uncertaintyStimulus
% slave PC
% yb, 06/15/2018, version 1.0

global Pstate

Pstate = struct; %clear it

Pstate.type = 'US';

Pstate.param{1} = {'predelay'  'float'      0       0                'sec'};
Pstate.param{2} = {'postdelay'  'float'     0       0                'sec'};
Pstate.param{3} = {'Nframes'  'int'     11       0                'frames'};

Pstate.param{4} = {'x_pos'       'int'      512       0                'pixels'};
Pstate.param{5} = {'y_pos'       'int'      410       0                'pixels'};
Pstate.param{6} = {'x_size'      'float'     15       1                'deg'};
Pstate.param{7} = {'y_size'      'float'     15       1                'deg'};

Pstate.param{8} = {'background'      'int'   128       0                'grayscales'};
Pstate.param{9} = {'contrast'    'float'     100       0                '%'};

Pstate.param{10} = {'stim_on'   'int'   200       0               'ms'};
Pstate.param{11} = {'stim_off'  'int'   100      0                'ms'};

Pstate.param{12} = {'image_folder'  'string'   '/Users/nauhauslab/Desktop/uncertainty_stimulus'       0                ''};

Pstate.param{13} = {'redgain' 'float'   1       0             ''};
Pstate.param{14} = {'greengain' 'float'   1       0             ''};
Pstate.param{15} = {'bluegain' 'float'   1       0             ''};

Pstate.param{16} = {'redbase' 'float'   .5       0             ''};
Pstate.param{17} = {'greenbase' 'float'   .5       0             ''};
Pstate.param{18} = {'bluebase' 'float'   .5       0             ''};

Pstate.param{19} = {'rseed'    'int'   1       0                ''};

Pstate.param{20} = {'eye_bit'    'int'   0       0                ''};
Pstate.param{21} = {'Leye_bit'    'int'   1       0                ''};
Pstate.param{22} = {'Reye_bit'    'int'   1       0                ''};



% Pstate.param{42} = {'light_bit'    'int'   0       0                ''};
% Pstate.param{43} = {'light_type'    'int'   0       0                ''};