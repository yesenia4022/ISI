function setacqinfo(trial)

%Get acquisition structure and set global

%loadAnalyzer needs to be run first 

global Analyzer ACQinfo

im = Load2phImage(1,1);

ACQinfo.numberOfFrames = length(Analyzer.syncInfo{trial}.acqSyncs);

ACQinfo.linesPerFrame = size(im,1);

ACQinfo.pixelsPerLine = size(im,2);