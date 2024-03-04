function setacqinfo_mat(trial)

%Get acquisition structure and set global

global ACQinfo twophDATADIR AUE

filepath = [twophDATADIR AUE ' ' sprintf('%03d',trial) '.mat'];

%load(filepath,'state')

load(filepath,'fileHeader')
state = parseHeaderNew(fileHeader);


ACQinfo = state.acq;