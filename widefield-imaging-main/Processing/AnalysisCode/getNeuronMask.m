function [nID gID] = getNeuronMask

global maskS

%%%Get rid of the glia:

masklabel1 = bwlabel(maskS.bwCell{1},4);
masklabel2 = bwlabel(maskS.bwCell{2},4);
neuronmask = maskS.bwCell{1};
celldom = unique(masklabel1);
k1 = 1; k2 = 1;
for i = 2:length(celldom)  %Don't use neuropil
    id = find(masklabel1 == celldom(i));
    if sum(masklabel2(id))==0  %if there is no overlap
        nID(k1) = i;
        k1 = k1+1;
    else
        gID(k2) = i;
        k2 = k2+2;
        neuronmask(id) = 0;
    end
end
%%%

maskS.neuronmask = neuronmask;
