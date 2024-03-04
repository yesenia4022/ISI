function prcoverlap = getcontroloverlap(xshift,yshift)

%percentage overlap of each neuron

global maskS

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);

%%%Figure out which control neurons overlapped too much with other neurons
mask_control = zeros(size(maskS.bwCell{1}));
mask = maskS.bwCell{1};

mask_control(yshift+1:end,xshift+1:end) = mask(1:end-yshift,1:end-xshift);
%mask_control(1:end-yshift,1:end-xshift) = mask(yshift+1:end,xshift+1:end);

masklabel_shift = bwlabel(mask_control,4);
shiftdom = unique(masklabel);
clear prcoverlap
for i = 2:length(shiftdom)
    id = find(masklabel_shift(:) == shiftdom(i));
    prcoverlap(i-1) = sum(mask(id))/length(id);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%