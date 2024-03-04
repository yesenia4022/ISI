function Grandposplots3

%Ian Nauhaus

global ACQinfo G_RChandles cellS maskS DM MK TC Analyzer

%%%%

DM = struct; MK = struct; TC = struct; %Reset these

%Make/store another structure related to the maskS, 'MK'
MK.Ncell = length(cellS.kernAll);  

MK.masklabel = bwlabel(maskS.bwCell{1},4);
MK.celldom = unique(MK.masklabel);
[MK.nID] = getNeuronMask;

for p = 1:MK.Ncell
    [idcelly idcellx] = find(MK.masklabel == MK.celldom(MK.nID(p)));
    MK.CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end

%Get/store the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string') ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with dtau spacing, and end at an estimate of kernDel(2)

%Get/store functional domains
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
load([logfileroot Analyzer.M.anim '\' expt],'domains')

DM.oridom = domains.oridom;
DM.colordom = domains.colordom;
DM.taudom = taudom;

DM.xdom = domains.xdom;
DM.xdom = DM.xdom-mean(DM.xdom); %center
degpersamp = getparam('x_size')/(length(DM.xdom)-1);  
DM.xdom = DM.xdom*degpersamp;  %Put units in degrees



%%
%getRFfromRandPos2

getRFfromRandPos_mouse

%%
plotRandPosMaps;

%%
%RCcolorplots

%%
pairWiseRandPosanalysis2

%%
%localNeighbors

