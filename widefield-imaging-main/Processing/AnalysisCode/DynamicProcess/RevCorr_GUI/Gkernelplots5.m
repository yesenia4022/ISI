function [maporidiff mapslopemag] = Gkernelplots5

%Ian Nauhaus

%5 was made for the mice

global ACQinfo G_RChandles cellS maskS DM MK

%%%%

DM = struct; MK = struct; %Reset these

%Make/store another structure related to the maskS, 'MK'
MK.Ncell = length(cellS.muTime);  

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
domains = getSeqInfo(trialdom);

if ~isempty(domains{1})
    id = 1;
else
    id = 2;
end
    
DM.oridom = domains{id}.oridom;
DM.sfdom = domains{id}.sfdom;
DM.phasedom = domains{id}.phasedom;
DM.colordom = domains{id}.colordom;
DM.taudom = taudom;


%%
getTCfromRevCorr5

%%
pairWiseRCanalysis2

%%
[maporidiff mapslopemag] = plotRCmaps;

%%
%RCcolorplots

%%
%localNeighbors
%localMapSlope  %Needs to be one or the other
%%


