function configPstate_geoCalibration

global Pstate

Pstate = struct; %clear it

Pstate.param{1} = {'Dot_Spacing'         'int'        300       0                'pixels'};
Pstate.param{2} = {'DotSize'       'int'      3       0                'pixels'};
Pstate.param{3} = {'calib_check'       'int'      0       0                'bit'};

