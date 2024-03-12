function Stimulator

%%

%Add controller folders to the path
root_controller = 'C:\stimulus_control\';
addpath([root_controller 'COM_acquisition'])
addpath([root_controller 'COM_display'])
addpath([root_controller 'Calibration'])
addpath([root_controller 'DisplayCode'])
addpath([root_controller 'GUIs'])
addpath([root_controller 'formula functions'])
addpath([root_controller 'onlineAnalysis'])
addpath([root_controller 'sync_inputs'])

%Remove analysis folder from the path
root_WFanalysis = 'C:\ISI acquisition and analysis\';
rmpath([root_WFanalysis 'Processing\AnalysisCode'])
rmpath([root_WFanalysis 'Processing\ISI_Processing'])
rmpath([root_WFanalysis 'Processing\ISIAnGUI'])
rmpath([root_WFanalysis 'Processing\ISIAnGUI\general'])
rmpath([root_WFanalysis 'Processing\AnalysisCode\ContrastResp'])
rmpath([root_WFanalysis 'Processing\AnalysisCode\DynamicProcess'])
rmpath([root_WFanalysis 'Processing\AnalysisCode\DynamicProcess\RevCorr_GUI'])
rmpath([root_WFanalysis 'Processing\offlineMovementCorrection'])


%Initialize stimulus parameter structures
configurePstate('PG')
configureMstate
configureLstate

%Host-Host communication
configDisplayCom    %stimulus computer (Psychtoolbox)
config2photonCom    %two-photon computer (Scanbox)


%NI USB input for ISI acquisition timing from frame grabber
%configSyncInput  

%configEyeShutter

%Open GUIs
MW
Looper 
paramSelect
