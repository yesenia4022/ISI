function configPstate_perGrater

%Ian Nauhaus

%periodic grater

global Pstate

Pstate = struct; %clear it

Pstate.type = 'PG';

id=1;
Pstate.param{id} = {'predelay'  'float'      2       0                'sec'};       id = id+1;
Pstate.param{id} = {'postdelay'  'float'     2       0                'sec'};       id = id+1;
Pstate.param{id} = {'stim_time'  'float'     1       0                'sec'};       id = id+1;

Pstate.param{id} = {'x_pos'       'int'      600       0                'pixels'};  id = id+1;
Pstate.param{id} = {'y_pos'       'int'      400       0                'pixels'};  id = id+1;
Pstate.param{id} = {'x_size'      'float'      3       1                'deg'};     id = id+1;
Pstate.param{id} = {'y_size'      'float'      3       1                'deg'};     id = id+1;
Pstate.param{id} = {'mask_type'   'string'   'none'       0                ''};     id = id+1;
Pstate.param{id} = {'mask_radius' 'float'      6       1                'deg'};     id = id+1;
Pstate.param{id} = {'x_zoom'      'int'   1       0                ''};             id = id+1;
Pstate.param{id} = {'y_zoom'      'int'   1       0                ''};             id = id+1;

Pstate.param{id} = {'altazimuth'   'string'   'none'       0                ''};    id = id+1;
Pstate.param{id} = {'tilt_alt'   'int'   0       0                'deg'};           id = id+1;
Pstate.param{id} = {'tilt_az'   'int'   0      0                'deg'};             id = id+1;
Pstate.param{id} = {'dx_perpbis'   'float'   0       0                'cm'};        id = id+1;
Pstate.param{id} = {'dy_perpbis'   'float'   0      0                'cm'};         id = id+1;

Pstate.param{id} = {'background'      'int'   128       0                ''};       id = id+1;
Pstate.param{id} = {'contrast'    'float'     100       0                '%'};      id = id+1;

Pstate.param{id} = {'ori'         'int'        0       0                'deg'};     id = id+1;

Pstate.param{id} = {'separable'   'int'     0       0                'bit'};        id = id+1;
Pstate.param{id} = {'st_profile'  'string'   'sin'       0                ''};      id = id+1;
Pstate.param{id} = {'st_phase'         'float'        180       0                'deg'}; id = id+1;

Pstate.param{id} = {'s_freq'      'float'      1      -1                 'cyc/deg'}; id = id+1;
Pstate.param{id} = {'s_profile'   'string'   'sin'       0                ''};          id = id+1;
Pstate.param{id} = {'s_duty'      'float'   0.5       0                ''};         id = id+1;
Pstate.param{id} = {'s_phase'      'float'   0.0       0                'deg'};     id = id+1;

Pstate.param{id} = {'t_profile'   'string'   'sin'       0                ''};      id = id+1;
Pstate.param{id} = {'t_duty'      'float'   0.5       0                ''};         id = id+1;
Pstate.param{id} = {'t_period'    'int'       20       0                'frames'};  id = id+1;
Pstate.param{id} = {'t_phase'      'float'   0.0       0                'deg'};     id = id+1;

Pstate.param{id} = {'noise_bit'      'int'   0       0                ''};          id = id+1;
Pstate.param{id} = {'noise_amp'      'float'   100       0                '%'};     id = id+1;
Pstate.param{id} = {'noise_width'    'int'   5       0                'deg'};       id = id+1;
Pstate.param{id} = {'noise_lifetime' 'float'   10       0             'frames'};    id = id+1;
Pstate.param{id} = {'noise_type' 'string'   'random'       0             ''};       id = id+1;

Pstate.param{id} = {'redgain' 'float'   1       0             ''};                  id = id+1;
Pstate.param{id} = {'greengain' 'float'   1       0             ''};                id = id+1;
Pstate.param{id} = {'bluegain' 'float'   1       0             ''};                 id = id+1;
Pstate.param{id} = {'redbase' 'float'   .5       0             ''};                 id = id+1;
Pstate.param{id} = {'greenbase' 'float'   .5       0             ''};               id = id+1;
Pstate.param{id} = {'bluebase' 'float'   .5       0             ''};                id = id+1;
Pstate.param{id} = {'colormod'    'int'   1       0                ''};             id = id+1;

Pstate.param{id} = {'mouse_bit'    'int'   0       0                ''};            id = id+1;

Pstate.param{id} = {'eye_bit'    'int'   1       0                ''};              id = id+1;
Pstate.param{id} = {'Leye_bit'    'int'   1       0                ''};             id = id+1;
Pstate.param{id} = {'Reye_bit'    'int'   1       0                ''};             id = id+1;

%Plaid variables
Pstate.param{id} = {'plaid_bit'    'int'        0       0             ''};          id = id+1;
Pstate.param{id} = {'plaid_function'  'string'   'add'       0                ''};     id = id+1;
Pstate.param{id} = {'contrast2'    'float'     10       0                '%'};      id = id+1;
Pstate.param{id} = {'ori2'         'int'        90       0                'deg'};   id = id+1;
Pstate.param{id} = {'st_phase2'         'float'        180       0                'deg'}; id = id+1;
Pstate.param{id} = {'st_profile2'  'string'   'sin'       0                ''};     id = id+1;
Pstate.param{id} = {'s_freq2'      'float'      1      -1                 'cyc/deg'}; id = id+1;
Pstate.param{id} = {'s_profile2'   'string'   'sin'       0                ''};     id = id+1;
Pstate.param{id} = {'s_duty2'      'float'   0.5       0                ''};        id = id+1;
Pstate.param{id} = {'s_phase2'         'float'        0       0                'deg'}; id = id+1;
Pstate.param{id} = {'t_profile2'   'string'   'sin'       0                ''};     id = id+1;
Pstate.param{id} = {'t_duty2'      'float'   0.5       0                ''};        id = id+1;
Pstate.param{id} = {'t_phase2'         'float'        0       0                'deg'}; id = id+1;
Pstate.param{id} = {'t_period2'         'int'        20       0                'frames'}; id = id+1;
Pstate.param{id} = {'altazimuth2'   'string'   'none'       0                ''};    id = id+1;
Pstate.param{id} = {'midpoint2' 'float'   0       0             ''};                  id = id+1;

Pstate.param{id} = {'harmonics'      'string'   '1'       0           'vector'};    id = id+1;
Pstate.param{id} = {'phaseShuff'      'int'   0       0                'bit'};      id = id+1;
Pstate.param{id} = {'speed_bit'  'int'   1      0                'bit'};  

