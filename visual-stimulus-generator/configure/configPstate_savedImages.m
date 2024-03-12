function configPstate_savedImages

global Pstate

Pstate = struct; %clear it

Pstate.param{1} = {'predelay'  'float'      2       0                'sec'};
Pstate.param{2} = {'postdelay'  'float'     2       0                'sec'};
Pstate.param{3} = {'stim_time'  'float'     1       0                'sec'};

Pstate.param{4} = {'x_pos'       'int'      600       0                'pixels'};
Pstate.param{5} = {'y_pos'       'int'      400       0                'pixels'};
Pstate.param{6} = {'x_size'      'float'      3       1                'deg'};
Pstate.param{7} = {'y_size'      'float'      3       1                'deg'};
Pstate.param{8} = {'mask_type'   'string'   'none'       0                ''};
Pstate.param{9} = {'mask_radius' 'float'      6       1                'deg'};
Pstate.param{10} = {'x_zoom'      'int'   1       0                ''};
Pstate.param{11} = {'y_zoom'      'int'   1       0                ''};
Pstate.param{12} = {'NxTile'      'int'   1       0                ''};
Pstate.param{13} = {'NyTile'      'int'   1       0                ''};

Pstate.param{14} = {'background'      'int'   128       0                ''};
Pstate.param{15} = {'contrast'    'float'     100       0                '%'};

Pstate.param{16} = {'ori'         'int'        0       0                'deg'};

Pstate.param{17} = {'h_per'      'int'   3       0                'frames'};
Pstate.param{18} = {'t_duty'    'float'     .1       0                ''};

Pstate.param{19} = {'redgain' 'float'   1       0             ''};
Pstate.param{20} = {'greengain' 'float'   1       0             ''};
Pstate.param{21} = {'bluegain' 'float'   1       0             ''};

Pstate.param{22} = {'redbase' 'float'   .5       0             ''};
Pstate.param{23} = {'greenbase' 'float'   .5       0             ''};
Pstate.param{24} = {'bluebase' 'float'   .5       0             ''};

Pstate.param{25} = {'gray_bit' 'int'   0       0             ''};

Pstate.param{26} = {'colormod' 'int'   1       0             ''};

Pstate.param{27} = {'rseed'    'int'   1       0                ''};

Pstate.param{28} = {'blankProb'    'float'   0       0                ''};

Pstate.param{29} = {'eye_bit'    'int'   0       0                ''};
Pstate.param{30} = {'Leye_bit'    'int'   1       0                ''};
Pstate.param{31} = {'Reye_bit'    'int'   1       0                ''};

Pstate.param{32} = {'image_folder'   'string'   '/Users/iannauhaus/Desktop/Textures'     0         ''};
Pstate.param{33} = {'fileID'   'string'   '*.png'     0         ''};
Pstate.param{34} = {'fileID_block2'   'string'   '*.png'     0         ''};
Pstate.param{35} = {'block_bit'   'int'   0     0         ''};
Pstate.param{36} = {'block_Nimage'   'int'   4     0         ''};

Pstate.param{37} = {'LPfilt' 'int'   0       0             'binary'};
Pstate.param{38} = {'LPcutoff' 'float'   8       0             'cyc/deg'};

Pstate.param{39} = {'distortbit' 'int'   0       0             ''};
Pstate.param{40} = {'truncXpercent' 'float'   100       0             'percent'};
Pstate.param{41} = {'truncYpercent' 'float'   100       0             'percent'};

Pstate.param{42} = {'repeatTrial_bit' 'int'   0       0             ''};
